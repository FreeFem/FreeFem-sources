#ifndef MY_LEX_HPP_
#define MY_LEX_HPP_
extern bool lexdebug;

class mylex  { 
  public:
  typedef const char * Key;
  typedef pair<int,aType> Value;
  struct Keyless : binary_function<Key,Key, bool>
   { bool operator()(const Key& x, const Key& y) const{ return strcmp(x,y)<0;} };  
  typedef   map<const char *,Value,Keyless> MapMotClef;
  typedef   map<const char *, deque<string> ,Keyless >  MapMacroDef;
  typedef   map<string, string  >  MapMacroParam;
//  typedef   map<const char *,Value,Keyless> MapMotType;
  typedef   MapMotClef::const_iterator const_iterator;
  typedef   MapMotClef::iterator iterator;
//  typedef   MapMotType::const_iterator const_iterType;
  public: 
  int linenumber,charnumber;

  private:
  int level;
  bool withmacropara;
  char buf[1024];
  int typetoken;
  struct xxxx { 
   istream * f;
   istream * nf;
   int l;
   int cas;
   
   const char * filename; 
   xxxx() : l(0), f(0) , filename(0),cas(0),nf(0) {}   
   void  open(mylex *lexx,const char * ff) ;
   void  readin(mylex *lexx,const string & s,int cc) ;
   void close() ;
   };
  friend struct mylex::xxxx;
  
  xxxx pilesource[100];
  istream & source() const {return  * pilesource[level].f;} 
  ostream & cout ;
//  ostream & cerr ;
  MapMotClef  MotClef;
  list<MapMacroDef> *listMacroDef;
  list<MapMacroParam> *listMacroParam;
  public:
   mylex(ostream & out):
    charnumber(0),
    linenumber(1),
    cout(out),
    level(-1),
    withmacropara(false),
    listMacroDef(new list<MapMacroDef>),
    listMacroParam(0) {
    listMacroDef->push_front(MapMacroDef());};
  string token() const;
  void print(ostream &f) const; 
  int scan();
  int lineno(){return linenumber;}
  char * YYText() { return buf;}
 void dump(ostream & f ) ;    
  void erreur(char * s) {
    cerr << " Erreur line number" <<linenumber << ": " << s << endl;
    throw(ErrorCompile("lex:",linenumber)); }

  bool InMotClef  (aType & t, int & r) const ;
  void  Add(Key k,int r,aType t);
  void Check(bool b,Key k,const char * s) {
    if(b) {cout <<"Add " << s << "  null: " << k << endl;
    CompileError();}}
  void Add(Key k,int i) ;//    {Check(!i,k,"mot clef");Add(k,i,0); }
  void Add(Key k,aType t);//   {Check(!t,k,"type");Add(k,TYPE,t); }
  void AddF(Key k,aType t);//  {Check(!t,k,"type");Add(k,FUNCTION,t); }
  void AddSpace(Key k,aType t);//  {Check(!t,k,"type");Add(k,FUNCTION,t); }
  const char * filename() const { 
    if ( level >=0 ) 
      return  pilesource[level].filename ? pilesource[level].filename : " -- in macro -- ";
    return "-- unkown --";}
    
  void input(const char *  filename) ;
  void input(const string &str,int nbparam) ;   
  bool close() ;
       
  private: 
    int basescan();  
    bool SetMacro(int &ret);
    bool CallMacro(int &ret);
    bool ExpandParam(int &ret);
    char * match(int i);
   void ErrorScan(char * s) {
      throw(ErrorCompile(s,lineno(),YYText() ) );}    
    
     
} ;
extern mylex *zzzfff;
#endif
