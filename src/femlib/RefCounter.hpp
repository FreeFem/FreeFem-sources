#ifndef REFCOUNTER_HPP
#define REFCOUNTER_HPP

#include "showverb.hpp" 

class RefCounter; 
class baseCountPointer;


// ruse pour utiliser le prive c RefCounter de 
//  pas de syntaxe pour des friends template 
class baseCountPointer { protected:
 void  add(const  RefCounter * c)  const; 
 void  destroyPtr(const RefCounter *  c)  const;
};

class RefCounter {

  mutable int count;
  protected:
  virtual ~RefCounter() {}
  RefCounter() : count(0) {}
  public:
  int destroy() const { 
   if(this) {throwassert(count>=0);
        if ( count--==0) {
	            SHOWVERB( cout << "True  destroy of " << this <<  endl);
             delete this;
             return true;} 
        else{ SHOWVERB(cout << " no destroy count=" << count+1 << " " << this <<  endl);
              return false;}}
   else return false;} 
   void increment() const {count++;} 
   void decrement() const {count--;}
 friend   class baseCountPointer;
// private:
  RefCounter(const RefCounter &) : count(0) {} 
  void operator=(const RefCounter &) { count=0;}
 

};

inline void baseCountPointer::add(const RefCounter * c)  const 
   { if (c) c->count++;}    
inline void baseCountPointer::destroyPtr(const RefCounter *   c)  const 
   { if (c) c->destroy();}    
 
template<class T>
class CountPointer: private baseCountPointer {
 T * c;
 public: 
 CountPointer() : c(0) {}
 CountPointer( T * a,bool master=false) :c(a) { if(!master) add(c);} 
 CountPointer(  T & a) :c(&a) { add(c);}
 CountPointer(const CountPointer & a) :c(a.c) { add(c);}
 ~CountPointer()  { destroyPtr(c);c=0;}
 //void destroy() const { destroyPtr(c);}
 void destroy()  { destroyPtr(c);c=0;}
 operator  T * ()   const { return c;}
 operator   T & () const  {return *c;}
  T& operator*() const {return *c;}
  T* operator->() const {return c;}
 bool operator==(const  CountPointer & n) const {return  n.c ==c;}
 bool operator!=(const  CountPointer & n) const {return  n.c !=c;}
 bool operator!() const { return !c;}
 void operator=(const  CountPointer & n) {
  if(*this != n) { destroyPtr(c); 
                c=n.c;
                add(c);               
              }}
                           
 void operator=( T * t) {
  if( c != t) { if(c) destroyPtr(c); 
                c=t;
                add(c);               
              }} 
//  for the compile time               
 void init() {c=0;} // 
 void master(T *t) {
     destroyPtr(c); 
     c=t;}
};

#endif
