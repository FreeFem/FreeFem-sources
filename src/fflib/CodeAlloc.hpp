
class CodeAlloc { public:

  static size_t nb,nbt,lg, nbdl,nbpx, chunk ;   
  static CodeAlloc ** mem;
  static bool cleanning;
  static void * lgmax;
  static bool sort;
  static bool isdel(int i)
  {
    return  ((char *) (void *)  mem[i] - (char *) 0) % 2 == 1;
  }
      
  static void setdel(int i)
  {
    mem[i] = (CodeAlloc *) (void *) (((char *) mem[i] - (char *) 0)+ 1);
  }
  static void resize(); 
  
 static  void * Add2CleanAtEnd(void * p)
  {
    if(p) {
      if(nbt>=nbpx) resize();
      if(nbt>0) sort = sort && mem[nbt-1] < p;
      nb++; 
      mem[nbt++]=(CodeAlloc*)p;  }
    return p;
  }
  
  void *operator new(size_t ll ) {
    lg+=ll;
    return Add2CleanAtEnd(::operator new(ll));} 

  
  static void Sort_mem();
  static  void clear();
  static void ErrorDel(void *pp);
  void operator delete(void * pp);
  virtual ~CodeAlloc() {}
  
};

template<class T> class CodeAllocT: public CodeAlloc{
  T * p;
  public: 
  CodeAllocT(int n): p(new T[n]) { assert(p);}
  static T * New(int n) { return (new CodeAllocT(n))->p;}
  ~CodeAllocT(){ delete [] p;}
};
