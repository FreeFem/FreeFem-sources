#ifndef STRING_HPP_
#define STRING_HPP_


#include <map>


template<class T>
string * toString(const T& a)
{ ostringstream r;
  r << a  ENDS ;
  return new string(r.str());
}


inline string * toStringCconst(const char * const &a)
{ return new string(a);
}
inline string * toStringC( char * const &a)
{ return new string(a);
}

template<class T>
string * PtoString(const  T * a)
{ ostringstream r;
  r << *a ENDS;
  return new string(r.str());
}

template<class T>
AnyType PtoStringA(void *, const AnyType &a)
{ 
  ostringstream r;
  r << *GetAny<T*>(a) ENDS ;
  return SetAny<string *>(new string(r.str()));
}
template<class T>
AnyType toStringA(void *, const AnyType &a)
{ 
  ostringstream r;
  r << GetAny<T>(a) ENDS ;
  return SetAny<string *>(new string(r.str()));
}

class String {  
  string  * p;
  public: 
//  String( string & pp) : p(&pp) {}
  String() : p(new string()) {/*cout << "String" << p <<","<<  *p << endl;*/}
//  String( string * c) : p(c) {cout << "String" << p <<","<< *p << endl;} 
  String(const String & c) : p(new string(c)) {/*cout << "String" << p <<","<< *p << endl;*/} 
  String(const string & c) : p(new string(c)) {/*cout << "String" << p <<","<< *p << endl;*/} 
  String(const char *  c) : p(new string(c)) {/*cout << "String" << p <<","<< *p << endl;*/} 
  String(const long & c) : p(toString(c)){/*cout << "String" << p <<","<< *p << endl;*/} 
  String(const double & c) : p(toString(c)){/*cout << "String" << p <<","<< *p << endl;*/} 
  String(const bool & c) : p(toString(c)){/*cout << "String" << p <<","<< *p << endl;*/} 
  String(const  long * c) : p(PtoString(c)){/*cout << "String" << p <<","<< *p << endl;*/} 
  String(const double * c) : p(PtoString(c)){/*cout << "String" << p <<","<< *p << endl;*/} 
  String & operator=(const String & s){delete p;p=new string(s);return *this;}
  String  operator+(const String & s)const {return String(new string(*p+*s.p));} 
  ~String(){/* cout << "~String" << p << *p << endl;*/delete p; p=0;}
   operator const string & () const {return *p;}
   operator  string & ()  {return *p;}
  friend inline  ostream & operator<<(ostream & f,const String & s) {throwassert(s.p); f << *s.p ; return f;}
  bool operator<(const String &t) const {assert(p && t.p);return *p<*t.p;} // correction FH feb 2004
  bool operator>(const String &t) const {assert(p && t.p);return *p>*t.p;} // correction FH feb 2004
};


template<class K,class V>
class MyMap {
  public:
  map<K,V> *m;
  
  MyMap() : m(new map<K,V>) {/*cout << "new MyMap:: " << m << endl;*/} 
  MyMap &operator=(MyMap &M){ 
  //   cout << " MyMap::= " << m << " = " << M.m << endl;
     delete m;m=new map<K,V>(M);
    }
  V &operator[](const K & k) {
  throwassert(m);
  typename map<K,V>::iterator i=m->find(k);
//  cout << k << " "  << " end? "  << (i==m->end())  << endl;
//  for(  map<K,V>::iterator ii=m->begin(); ii != m->end();ii++)
 //    cout << " MyMap ::   m="<< m << ": " <<  ii->first  << "  -> " <<  ii->second <<endl;
 // map<K,V>::iterator j=m->find(k);
  
 // cout << " m->find(k)->second   " <<  i->second  << ";" <<  j->second <<endl;
   
  if  (i==m->end())
     i=m->insert(map<K,V>::value_type(k,V())).first;
  V &  v= i->second ;  
  return v;
}
  ~MyMap(){delete m;/*cout << "MyMap:: delete "<< m << endl;*/m=0;} 
  private:
    MyMap(const MyMap &M):m(new map<K,V>(*M.m)) {}
 
};
template<class A,class B>
struct pairless : binary_function<pair<A,B>,const char *, bool>
   { 
    typedef pair<A,B> Key;
    bool operator()(const Key& x, const Key& y) const { return  x.first<y.first ? true
       :  (  (x.first == y.first)  ? (x.second < y.second) : false );} };
       
typedef   map<pair<aType,aType>,aType,pairless<aType,aType> > Map_type_of_map; 
    
extern  Map_type_of_map map_type_of_map ; //  to store te type 
//  of a map  of  A[B] 


#endif
