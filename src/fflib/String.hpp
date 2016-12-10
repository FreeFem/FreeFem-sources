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
#ifndef STRING_HPP_
#define STRING_HPP_


#include <map>
#include <cstdio>

// BUG option compilation -fast 
/*
template<class T>
string * toString(const T& a)
{ ostringstream r;
  r << a  ENDS ;
  return new string(r.str());
}

*/
//  to be sure new and delele be in see dll for windows
 string  *newstring();
 string  *newstring(const string & c);
 string  *newstring(const char * c);
 void freestring(const string * c);
//fh June 2016 ....  Hard bug 
inline string * toString(const double& a)
{
  char buf[30];
  sprintf(buf,"%g",a);
 return newstring(buf);
}
inline string * toString(const long& a)
{
  char buf[30];
  sprintf(buf,"%ld",a);
  return newstring(buf);
}
inline string * toString(const bool& a)
{
  return newstring(a?"T":"F");
}
inline string * toString(const complex<double> & a)
{
  char buf[60];
  sprintf(buf,"%g%+gi",a.real(),a.imag());
  return newstring(buf);
}


inline string * toStringCconst(const char * const &a)
{ return newstring(a);
}
inline string * toStringC( char * const &a)
{ return newstring(a);
}

template<class T>
string * PtoString(const  T * a)
{ ostringstream r;
  r << *a ENDS;
  return newstring(r.str());
}

template<class T>
AnyType PtoStringA(void *, const AnyType &a)
{ 
  ostringstream r;
  r << *GetAny<T*>(a) ENDS ;
  return SetAny<string *>(newstring(r.str()));
}
template<class T>
AnyType toStringA(void *, const AnyType &a)
{ 
  ostringstream r;
  r << GetAny<T>(a) ENDS ;
  return SetAny<string *>(newstring(r.str()));
}

class String {
  
  string  * p;
  public:

//  String( string & pp) : p(&pp) {}
  String() : p(newstring()) {/*cout << "String" << p <<","<<  *p << endl;*/}
  void init() { p= newstring();} // Add FH march 2010
void destroy() { freestring(p);p=0;} // Add FH march 2010
//  String( string * c) : p(c) {cout << "String" << p <<","<< *p << endl;} 
  String(const String & c) : p(newstring(*c.p)) {/*cout << "String" << p <<","<< *p << endl;*/}
  String(const string & c) : p(newstring(c)) {/*cout << "String" << p <<","<< *p << endl;*/}
  String(const string * c) : p(newstring(*c)) {/*cout << "String" << p <<","<< *p << endl;*/}
  String(const char *  c) : p(newstring(c)) {/*cout << "String" << p <<","<< *p << endl;*/}
  String(const long & c) : p(toString(c)){/*cout << "String" << p <<","<< *p << endl;*/} 
  String(const double & c) : p(toString(c)){/*cout << "String" << p <<","<< *p << endl;*/} 
  String(const bool & c) : p(toString(c)){/*cout << "String" << p <<","<< *p << endl;*/} 
  String(const  long * c) : p(PtoString(c)){/*cout << "String" << p <<","<< *p << endl;*/} 
  String(const double * c) : p(PtoString(c)){/*cout << "String" << p <<","<< *p << endl;*/} 
  String & operator=(const String & s){ freestring(p);p=newstring(s);return *this;}
  String  operator+(const String & s)const {return String(newstring(*p+*s.p));} 
  ~String(){ if(verbosity>999999) cout << "~String" << p <<" "<<  *p << "." << endl; freestring(p); p=0;}
   operator const string & () const {return *p;}
   operator  string & ()  {return *p;}
   operator  const string * ()  const {return p;}
   operator  string *& ()   {return p;}
    
   string **  getap()  {return &p;}
  friend inline  ostream & operator<<(ostream & f,const String & s) {throwassert(s.p); f << *s.p ; return f;}
  bool operator<(const String &t) const {assert(p && t.p);return *p<*t.p;} // correction FH feb 2004
  bool operator>(const String &t) const {assert(p && t.p);return *p>*t.p;} // correction FH feb 2004
};


template<class K,class V>
class MyMap {
  public:
    typedef V Value;
    typedef K Key;
    typedef map<K,V> MAP;
    typedef typename map<K,V>::iterator iterator;
  map<K,V> *m;
  
  MyMap() : m(new map<K,V>) {/*cout << "new MyMap:: " << m << endl;*/} 
  MyMap &operator=(MyMap &M){ 
  //   cout << " MyMap::= " << m << " = " << M.m << endl;
     delete m;m=new map<K,V>(M);
    }
    bool exist(const K & k) const { return m&&  (m->find(k)!= m->end());}
   bool  insert(const K & k,const V & v)
    {
        typedef typename map<K,V>::value_type   value_type;
        return m->insert(value_type(k,v)).second;
    }
    bool insert(const K & k, V * const v)
    {
        typedef typename map<K,V>::value_type   value_type;
        return m->insert(value_type(k,*v)).second;
    }
    
  V &operator[](const K & k) {
  throwassert(m);
  typename map<K,V>::iterator i=m->find(k);
//  cout << k << " "  << " end? "  << (i==m->end())  << endl;
//  for(  map<K,V>::iterator ii=m->begin(); ii != m->end();ii++)
 //    cout << " MyMap ::   m="<< m << ": " <<  ii->first  << "  -> " <<  ii->second <<endl;
 // map<K,V>::iterator j=m->find(k);
  
 // cout << " m->find(k)->second   " <<  i->second  << ";" <<  j->second <<endl;
  typedef typename map<K,V>::value_type   value_type;
  if  (i==m->end())
     i=m->insert(value_type(k,V())).first;
  V &  v= i->second ;  
  return v;
}
    ~MyMap(){if(verbosity>99999)cout << " ~MyMap:: delete "<< m << endl; delete m;m=0;}
    void destroy(){if(verbosity>99999) cout << " MyMap:: destroy "<< m << endl;delete m;m=0;}
    void init() {  m=new MAP; if(verbosity>99999) cout << " MyMap:: int  "<< m << endl; }
  private:
    MyMap(const MyMap &M):m(new map<K,V>(*M.m)) {}
public:
    friend inline  ostream & operator<<(ostream & f,const MyMap & s)
    {
        if(s.m)
        for(MyMap::iterator i=s.m->begin(); i!= s.m->end(); ++i)
            f << i->first << " " << i->second << endl;
        return f;}

    
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
