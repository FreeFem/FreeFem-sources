#include <queue>
// class for 1 linear form
// just copy a array 

//  Warning this class are use at compilating time
//  ----------------------------------------------
using Fem2D::operatortype;
using Fem2D::op_id;
using Fem2D::op_dx;
using Fem2D::op_dy;
using Fem2D::last_operatortype;

template<class T> 
inline  T * NewCopy(const  T * const o,int n)
{ 
  int  m=n;
  T * c = new T [n];
  for (int i=0;i<m;i++)
    c[i]=o[i];
  return c;
}
template<class T> 
inline  T * NewCopy(const  T * const o,int n,int m)
{ 
  throwassert(m<=n);
  T * c = new T [n];
  for (int i=0;i<m;i++) c[i]=o[i];
  return c;
}

template <class T1, class T2, class T3>
struct triplet
{
	typedef T1 first_type;
	typedef T2 second_type;
	typedef T3 third_type;
	T1 first;
	T2 second;
	T3 third;
	triplet() : first(),second(),third()   {}
	triplet(const T1& x, const T2& y, const T3& z)  : first(x),second(y), third(z) {}
	template<class U, class V,class W> inline
	triplet(const triplet<U, V, W>& p) : first(p.first),second(p.second),third(p.third) {}
};

template <class T1, class T2,class T3>
inline triplet<T1, T2, T3> make_triplet(T1 x, T2 y, T3 z)
{
	return triplet<T1, T2, T3>(x, y, z);
}


template<class I,class R>
 class LinearComb : public E_F0mps { public: 
   typedef I TI;
   typedef R TR;
   typedef size_t size_type;
   typedef pair<I,R> K;
   typedef vector<K> array;
   typedef typename array::const_iterator const_iterator;
   typedef typename array::iterator iterator;
   array v;
   vector<size_type> where_in_stack_opt;
   Expression optiexp0,optiexpK;
   
   LinearComb(): v(),optiexp0(),optiexpK(),where_in_stack_opt(0) {}
   LinearComb(const I& i,const R& r) :v(1),optiexp0(),optiexpK(),where_in_stack_opt(0) {
    v[0]=make_pair<I,R>(i,r);}
    
   LinearComb(const LinearComb &l) 
      :v(l.v),optiexp0(l.optiexp0),optiexpK(l.optiexpK),where_in_stack_opt(l.where_in_stack_opt){}  
       
   void operator=(const LinearComb<I,R> &l) {v=l.v;}
   
   const I * simple() const { if (v.size()==1) return & v.begin()->first;else return  0;}     
   void  add(const I& i,const R &r)  { 
     for (iterator k=v.begin();k!=v.end();k++)
       if (k->first == i) {k->second += r;return ;}
     v.push_back(make_pair<I,R>(i,r));     
   }
   
   size_type size() const { return v.size();}   
   
   const K & operator[](size_type i) const { return v[i];}
   
   void operator+=(const LinearComb & l) {
     for (const_iterator k=l.v.begin();k!=l.v.end();k++)
       { const K & kk(*k);
       add(kk.first,kk.second);}
   }
   
   void operator*=(const R & r) {
     for (iterator k=v.begin();k!=v.end();k++)
       {K & kk(*k);
       kk.second = kk.second*r;}
   }     
   
   void operator/=(const R & r) {
     for (iterator k=v.begin();k!=v.end();k++)
       {K & kk(*k);
       kk.second = kk.second/r; }
   }     
   
  AnyType operator()(Stack )  const {
    return SetAny<const LinearComb<I,R> * >(this);}
  operator aType () const { return atype<const LinearComb >();}         
    
    
  bool mappable(bool (*f)(const  R &)) const {
     for (const_iterator k=v.begin();k!=v.end();k++)
       if (!(*f)(k->second)) return false;
      return true;}
      
  void mapping(R (*f)(const R &)) {
     for (iterator k=v.begin();k!=v.end();k++)
       k->second=(*f)(k->second) ;}
  
  int MaxOp() const {
    int m=0;
     for (const_iterator k=v.begin();k!=v.end();k++)
       m= maxop(m,(k->first)) ;
     return m;}
    
 void DiffOp(KN_<bool> &d) const {
     assert(d.N() >= last_operatortype);
     d=false;
     for (const_iterator k=v.begin();k!=v.end();k++)
       SetOp(d,k->first);
       
   
  }  
  
  ostream & dump(ostream &f) const {
    int n = size();
    for (int i=0; i<n; i++)
      {
        const K & ri=v[i];
        Expression ee= ri.second.LeftValue();
        f  << "\n\t\t"<< i << " " << ri.first << ": "  ;
       f << " :  type exp: " << typeid(*ee).name() << " "<<endl;
      }

    return f;
  } 
 LinearComb * Optimize(Block * b) const 
  {
    const bool kdump=false;
    if (kdump)
    cout << "\n\n Optimize " << endl;
    LinearComb * r=new LinearComb(*this);
    LinearComb  &rr= *r;
    int n = rr.size();
    
    deque<pair<Expression,int> > ll;
    MapOfE_F0 m;
    rr.where_in_stack_opt.resize(n);
    size_type top = b->OffSet(0), topbb=top; // FH. bofbof ??? 
    for (int i=0; i<n; i++)
      {
      const K & ri=rr.v[i];
       Expression ee= ri.second.LeftValue();
       if (kdump)
       cout << "Optimize :  type exp: " << typeid(*ee).name() << " "<<endl;
       rr.where_in_stack_opt[i]=ee->Optimize(ll, m, top);
       if (kdump)
       cout  << "\n\t\t"<< i << " " << ri.first << ": " << rr.where_in_stack_opt[i] << endl;
      }
      
    b->OffSet(top-topbb);
    //  
    int k=ll.size(),k0=0,k1=0;
    for (int i=0;i<k;i++)
        if (ll[i].first->MeshIndependent()) k0++;
    deque<pair<Expression,int> > l0(k0),l1(k-k0);
    k0=0,k1=0;
    for (int i=0;i<k;i++)
       if (ll[i].first->MeshIndependent()) 
         {
          if (kdump)
          cout << " mi " << ll[i].second << " " << *(ll[i].first) << endl;
          l0[k0++]=ll[i];
         }
        else 
         {
          if (kdump)
          cout << " md " << ll[i].second << " " << *(ll[i].first) << endl;
          l1[k1++]=ll[i];
         }
    if (k0)      
      rr.optiexp0 = new E_F0_Optimize(l0,m,0);  
    if (k1) 
      rr.optiexpK = new E_F0_Optimize(l1,m,0);
    if (kdump) cout << "LinearCom Optimize k0(mi) = " << k0 << " k1 = " << k1 << "\n\n"<<endl;
    return r;
  }    
   
};



template<class I,class R>
LinearComb<I,R> operator+(const LinearComb<I,R> & a,const LinearComb<I,R> & b)
 {LinearComb<I,R> r(a);r+=b;return r;}

template<class I,class R>
LinearComb<I,R> operator*(const LinearComb<I,R> & a,const R & b)
 {LinearComb<I,R> r(a);r*=b;return r;}
template<class I,class R>
LinearComb<I,R> operator*(const R & b,const LinearComb<I,R> & a)
 {LinearComb<I,R> r(a);r*=b;return r;}
 

 
class MGauche :public pair<int,operatortype> {public:
  MGauche() {}
  MGauche(int i,operatortype j) {first = i;second= j;}
  MGauche(const pair<int,operatortype> &p) : pair<int,operatortype>(p){}
  bool operator==(const MGauche& a) const {
    return static_cast<bool>(first == a.first && second == a.second);}
  int maxop(int op) const { return Max(op,(int) second);}
    
};

class MDroit :public pair<int,operatortype> {public:
  MDroit(){}
  MDroit(int i,operatortype j) {first = i;second =j;}
  //   first : number of unknow 
  // second  : number of operator ( 0 Id, 1 dx, 2 dy)
  MDroit(const pair<int,operatortype> &p) : pair<int,operatortype>(p){}
  bool operator==(const MDroit& a) const {
    return static_cast<bool>(first == a.first && second == a.second);}
};

inline ostream & operator<<(ostream & f,const MDroit & p)
{ f << p.first <<','<<p.second ;
  return f;}
inline ostream & operator<<(ostream & f,const MGauche & p)
{ f << p.first <<';'<<p.second ;
  return f;}
  
template<typename A,typename B>
inline ostream & operator<<(ostream & f,const pair<A,B> &p)
{ f << p.first <<" "<<p.second ;
  return f;}
  
  
//extern const C_F0 & One, &Zero;
C_F0 operator*(const C_F0 &,const C_F0 &);
C_F0 & operator+=(const C_F0 &,const C_F0 &);
 
typedef LinearComb<MGauche,C_F0> LinearOperatorG;
typedef LinearComb<MDroit ,C_F0> LinearOperatorD;

typedef LinearComb<pair<MGauche,MDroit>,C_F0> BilinearOperator;

inline  int maxop(int op,const MGauche & v) 
    { return Max(op,(int) v.second);}
inline  int maxop(int op,const MDroit & v) 
    { return Max(op,(int) v.second);}
inline  int maxop(int op,const pair<MGauche,MDroit> & b) 
    { return Max(op,(int) b.first.second,(int) b.second.second );}
 
inline BilinearOperator operator*(const LinearOperatorG & a,const LinearOperatorD & b) 
 {
   BilinearOperator r;
   for (LinearOperatorG::const_iterator i=a.v.begin();i!=a.v.end();i++)
     for (LinearOperatorD::const_iterator j=b.v.begin();j!=b.v.end();j++)
       {
	 const LinearOperatorG::K  vi(*i);
	 const LinearOperatorD::K  vj(*j);
	
	 r.add(make_pair(vi.first,vj.first),vi.second*vj.second);
	   
       }
    return r;
 } 
inline BilinearOperator operator*(const LinearOperatorD & b,const LinearOperatorG & a) 
 {
   BilinearOperator r;
   for (LinearOperatorG::const_iterator i=a.v.begin();i!=a.v.end();i++)
     for (LinearOperatorD::const_iterator j=b.v.begin();j!=b.v.end();j++)
       {
	 const LinearOperatorG::K  vi(*i);
	 const LinearOperatorD::K  vj(*j);
	
	 r.add(make_pair(vi.first,vj.first),vi.second*vj.second);
	   
       }
    return r;
 } 
ostream & operator<<(ostream & f,const  BilinearOperator & a); 
 
inline ostream & operator<<(ostream & f,const  BilinearOperator & a)
{
  int k=0;
  for (BilinearOperator::const_iterator i=a.v.begin();i!=a.v.end();i++)
    {
      const BilinearOperator::K  vi(*i);
      char * www[]={" ","_x ","_y "};
      const  pair<pair<int,int>,pair<int,int> > i1(vi.first);
      const  pair<int,int> ii(i1.first),jj(i1.second);
      f << *(const E_F0 *) vi.second <<  char('u'+ii.first) << www[ii.second] << " " << char('u'+jj.first)<<"'" << www[jj.second] ;
      if ( (++k%5)==0) f << endl ; else f << " ";
    }
  return f;
}
typedef LinearOperatorD LOperaD;
typedef LinearOperatorG LOperaG;
typedef BilinearOperator Opera;
inline LOperaG DiffG(int i,operatortype j) { return LOperaG(make_pair(i,j),*pOne);}
inline LOperaD DiffD(int i,operatortype j) { return LOperaD(make_pair(i,j),*pOne);}

inline LOperaG *newU_(int i) {  LOperaG * r;
  r= new LOperaG(make_pair(i,op_id),*pOne);
 SHOWVERB( cout << "newU_ " << r << endl);
  return r; }
inline LOperaG *newU_x(int i) { return new LOperaG(make_pair(i,op_dx),*pOne);}
inline LOperaG *newU_y(int i) { return new LOperaG(make_pair(i,op_dy),*pOne);}
inline LOperaD *newV_(int i) { return new LOperaD(make_pair(i,op_id),*pOne);}
inline LOperaD *newV_x(int i) { return new LOperaD(make_pair(i,op_dx),*pOne);}
inline LOperaD *newV_y(int i) { return new LOperaD(make_pair(i,op_dy),*pOne);}

template<class L>
 L *  Diff(const L *   u,const  operatortype & d) { 
     throwassert(u);
     L * r= new  L(*u);
    for (typename L::iterator i=r->v.begin();i!=r->v.end();i++)
      {   operatortype & dd=i->first.second;
          
         if (dd != op_id)
           { throwassert(0); i->first.second = op_id;  } // a faire          
         else {
           throwassert(i->second.EvaluableWithOutStack());// a faire derivation des fonctions
           dd = d ; }
    }
    return r;}
    


template<class L>
int  MaxOp(const L *   u) { 
     throwassert(u);
     int op=0;
    for (typename L::const_iterator i=u->v.begin();i!=u->v.end();i++)          
          op = maxop(op,i->first);
    return op;}
    
    
    
template<class L,class A,L* ff(const basicAC_F0 & args)>
class CODE_L1 { public:
  typedef L* Result;
  typedef L* (*func)(const basicAC_F0 & args) ;
  static L* f(const basicAC_F0 & args) { return ff(args);} //  ff :A-> L*
  static ArrayOfaType  typeargs() {return ArrayOfaType(atype<A>);}
};

template<class L>
class CODE_L_Add { public:
  typedef const L* Result;
  static  E_F0 * f(const basicAC_F0 & args) { 
    const L * a(dynamic_cast<const L*>((Expression) args[0]));
    const L * b(dynamic_cast<const L*>((Expression) args[1]));
    throwassert(a  && b );
    return new L(*a+*b);}
   static ArrayOfaType  typeargs() {return ArrayOfaType(atype<const L*>(),atype<const L*>());}
};
template<class L>
class CODE_L_Sub { public:
  typedef const L* Result;
  static  E_F0 * f(const basicAC_F0 & args) { 
    const L * a(dynamic_cast<const L*>((Expression) args[0]));
    const L * b(dynamic_cast<const L*>((Expression) args[1]));
    throwassert(a  && b );
    L * bb = new L(*pminusOne * *b);
    return new L(*a+*bb);}
   static ArrayOfaType  typeargs() {return ArrayOfaType(atype<const L*>(),atype<const L*>());}
};

template<class L>
class CODE_L_Minus { public:
  typedef const L* Result;
  static  E_F0 * f(const basicAC_F0 & args) { 
    const L * a(dynamic_cast<const L*>((Expression) args[0]));
    throwassert(a  && pminusOne );
    return new L(*pminusOne * *a);}
   static ArrayOfaType  typeargs() {return ArrayOfaType(atype<const L*>());}
};

template<class L,class A,class B>
class CODE_L_Mul { public:
  typedef const L* Result;
  static  E_F0 * f(const basicAC_F0 & args) { 
    const A * a(dynamic_cast<const A*>((Expression) args[0]));
    const B * b(dynamic_cast<const B*>((Expression) args[1]));
    SHOWVERB(cout << " CODE_L_Mul " << a << " " << b << endl);
    throwassert(a  && b );
    return new L(*a * *b);}
   static ArrayOfaType  typeargs() {return ArrayOfaType(atype<const A*>(),atype<const B*>());}
};

template<class T,class L>
class CODE_L_MulRL { public:
  typedef const L* Result;
  static  E_F0 * f(const basicAC_F0 & args) { 
    const L * b(dynamic_cast<const L*>((Expression) args[1]));
    if (!b) {
      
      cout << " --- " << ((Expression) args[1]) << typeid((Expression) args[1]).name() << endl;
    throwassert( b );
    }
    return new L(to<T>(args[0]) * *b);}
   static ArrayOfaType  typeargs() {return ArrayOfaType(atype<T>(),atype<const L*>());}
};

template<class L,class T>
class CODE_L_MulLR { public:
  typedef const L* Result;
  static  E_F0 * f(const basicAC_F0 & args) { 
    const L * a(dynamic_cast<const L*>((Expression) args[0]));
    throwassert( a );
    return new L(to<T>(args[1]) * *a);}
   static ArrayOfaType  typeargs() {return ArrayOfaType(atype<Result>(),atype<T>());}
};

template<class L,class T>
class CODE_L_DivLR { public:
  typedef const L* Result;
  static  E_F0 * f(const basicAC_F0 & args) { 
    const L * a(dynamic_cast<const L*>((Expression) args[0]));
    throwassert( a );
    return new L(C_F0(TheOperators,"/",*pOne,args[1]) * *a);}
   static ArrayOfaType  typeargs() {return ArrayOfaType(atype<Result>(),atype<T>());}
};

template<class L,operatortype op>
class CODE_Diff { public:
  typedef const L* Result;
  static  E_F0 * f(const basicAC_F0 & args) { 
   const L * a=dynamic_cast<const L*>((Expression) args[0]) ;
    throwassert( a );
    return  Diff(a,op);}
   static ArrayOfaType  typeargs() {return ArrayOfaType(atype<Result>());}
};

enum TheCode_VF { Code_Jump=1, Code_Mean=2, Code_OtherSide=3};

template<class L,TheCode_VF code>
class Code_VF { public:  
  typedef const L* Result;
  static  E_F0 * f(const basicAC_F0 & args) { 
   const L * u=dynamic_cast<const L*>((Expression) args[0]) ;
    assert( u );
     L * r= new  L(*u);
    for (typename L::iterator i=r->v.begin();i!=r->v.end();i++)
      {   operatortype & dd=i->first.second;
          assert(dd<last_operatortype);
          dd = (operatortype) ((int) dd+last_operatortype*code) ; 
     }    
    return r;}
   static ArrayOfaType  typeargs() {return ArrayOfaType(atype<Result>());}
};
         
