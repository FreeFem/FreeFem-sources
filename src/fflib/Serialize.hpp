#ifndef SERIALEZE_HPP_
#define SERIALEZE_HPP_

struct MPIrank;
 class Serialize {  
   // we store a refcounter in the pointer p a adresse p-sizeof(long)
   // so we can use the copy constructor
   private: 
   
  size_t lg;
  const char *what;
  char * p; 
  public: 
  Serialize(size_t lgg,const char * wht): 
    lg(lgg), what(wht) , p((new char[lg+sizeof(long)])+sizeof(long)) 
   { cout << " begin count()=0 " << endl;
    count()=0; }
  
  ~Serialize(){ if(count()--==0) delete [] (p-sizeof(long));}
  size_t size() const { return lg;}
  //  mpi routine
  void mpisend(const MPIrank &,long tag);
  Serialize(const MPIrank &,const char * wht,long tag); 
  //  end mpi routine 
  operator void *() { return p;} 
  operator char *() { return p;} 
  bool samewhat(const char * w) const { return strncmp(what,w,8)==0; }
  Serialize(const Serialize & s) :lg(s.lg),p(s.p),what(s.what) { count()++; } 
  
 template<typename T>  inline void get(size_t & k,T & x) const
   { 
     assert(k<=lg+sizeof(T));
     memcpy(&x,p+ k,sizeof(T));
     k +=  sizeof(T);   }
 template<typename T>  inline void put(size_t & k,const T & x) 
   { 
     if ( !(k<=lg+sizeof(T)) )
     {
       cout << " assert put " << k << " <=" << lg + sizeof(T) << endl;
       assert((k<=lg+sizeof(T)));
     }
    memcpy( p + k,&x,sizeof(T));
   k += sizeof(T);  }
   
  private:
   long & count() const  { return * (long*) (void*) (p-sizeof(long));}
   void operator=(Serialize & s) ; // no affectation

 };
#endif
