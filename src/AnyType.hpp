//   definition of the a type to store any bacis type
//   void *, long, double, complex, all bigger type 
//   are pointer 
//#define WITHCHECK
class basicForEachType;
typedef const  basicForEachType * aType;

typedef  unsigned char  AnyData[24]; 

extern map<const string,basicForEachType *> map_type;

template <class T,bool ok>
  class CheckSize { } ;
template <class T>
class CheckSize<T,false> {
  CheckSize() { } // private constructor to see the error
} ;
class AnyTypeWithCheck  { 

 private:  
    union { 
      AnyData data;
      void * p;
      double r;  
      long l;
      bool b;
//      complex<double> c;
    };
 public:
    const basicForEachType *ktype; 
//    friend ostream & operator<<(ostream &f,const AnyTypeWithCheck &p);  
    AnyTypeWithCheck(): ktype(0){}
    AnyTypeWithCheck(long ll)   {l=ll;ktype=map_type[typeid(long).name()];}
    AnyTypeWithCheck(double dd) {r=dd;ktype=map_type[typeid(double).name()];}
    AnyTypeWithCheck(bool bb ) {b=bb;ktype=map_type[typeid(bool).name()];}
    
    void operator =(void *pp){p=pp;}
};

class AnyTypeWithOutCheck  { 
    
 public: 
    union { 
      AnyData data;
      void * p;
      double r;  
      long l;
      bool b;
//      complex<double> c;
    };
 public:
    AnyTypeWithOutCheck(){}
    AnyTypeWithOutCheck(long ll)   {l=ll;}
    AnyTypeWithOutCheck(double dd) {r=dd;}
    AnyTypeWithOutCheck(void *pp ) {p=pp;}
    AnyTypeWithOutCheck(bool bb ) {b=bb;}
  
    void operator =(void *pp){p=pp;}
};

#ifdef WITHCHECK
typedef   AnyTypeWithCheck AnyType; 
#else
typedef   AnyTypeWithOutCheck AnyType; 
#endif
static AnyType Nothing;

  
#ifdef WITHCHECK  

template<typename T>   AnyTypeWithCheck SetAny(const T & x)
  { 
	AnyTypeWithCheck any;
    CheckSize<T,sizeof(T)<= sizeof(AnyData) >();
	throwassert( (any.ktype=map_type[typeid(T).name()]));
	throwassert (sizeof(T) <= sizeof(AnyData));
	//a.t=x;
	memcpy(&any,&x,sizeof(x));
	return any;
  }
  
inline AnyTypeWithCheck PtrtoAny(void * p,aType r)
  {
	AnyTypeWithCheck any;
    any=p;
	throwassert(any.ktype=r);    
	return any;
  }

  
#else

 template<typename T>   AnyTypeWithOutCheck inline SetAny(const T & x)
  { 
   AnyTypeWithOutCheck any;
   CheckSize<T,sizeof(T)<= sizeof(AnyData) >();

  // memcpy(&any,&x,sizeof(x));
   any = *(  (AnyTypeWithOutCheck *) (void *) &x); 
//  toto.x=x;
   return any;
  }
  
inline AnyTypeWithOutCheck PtrtoAny(void * p,aType )
  {
    return p;
  }
 template<>  inline AnyTypeWithOutCheck SetAny<double>(const double & x)
  { return x;}
 template<> inline  AnyTypeWithOutCheck SetAny<long>(const long & x)
  { return x;}
 template<> inline  AnyTypeWithOutCheck SetAny<bool>(const bool & x)
  { return x;}
 template<typename T>  inline  AnyTypeWithOutCheck SetAny( T * const & x)
 { return (void *) x;}
  
  
#endif
 template<class T> inline const T& GetAny(const AnyTypeWithCheck & x);
 
 template<class T> inline  const T& GetAny(const AnyTypeWithCheck & x)
{ 
 CheckSize<T,sizeof(T)<= sizeof(AnyData) >();
 if (x.ktype!=map_type[typeid(T).name()])
   { cerr<< "GetAny: PB type <"<<typeid(T).name() << "> <=" << *x.ktype << endl;
   throw(ErrorExec("exit",1));}
 return *static_cast<const T*>(static_cast<const void*>(&x));
}





 template<typename T> inline  const T& GetAny(const AnyTypeWithOutCheck & x)
{ 
 CheckSize<T,sizeof(T)<= sizeof(AnyData) >();
 return *static_cast<const T*>(static_cast<const void*>(&x));
}
 template<> inline const double& GetAny<double>(const AnyTypeWithOutCheck & x)
 { return x.r;} 
 template<>  inline const long& GetAny<long>(const AnyTypeWithOutCheck & x)
 { return x.l;} 
 template<>  inline const bool& GetAny<bool>(const AnyTypeWithOutCheck & x)
 { return x.b;} 
 
//template<> template<class T>  inline   T * const & GetAny<T*>(const AnyTypeWithOutCheck & x)
// { return x.p;}


