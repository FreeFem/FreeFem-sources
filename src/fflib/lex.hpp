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
#ifndef MY_LEX_HPP_
#define MY_LEX_HPP_
// with New version of macro expansion more simple and more stable 
// FH jan 2005
#include <stack> 
#include "environment.hpp"
extern bool lexdebug;
extern long mpisize,mpirank;

class mylex : public CodeAlloc { 
  public:
  typedef const char * Key;
  typedef pair<int,aType> Value;
  struct Keyless : binary_function<Key,Key, bool>
   { bool operator()(const Key& x, const Key& y) const{ return strcmp(x,y)<0;} };  
  typedef   map<const char *,Value,Keyless> MapMotClef;
  typedef   map<const char *, deque<string> ,Keyless >  MapMacroDef;
  typedef   map<string, string  >  MapMacroParam;
  typedef   MapMotClef::const_iterator const_iterator;
  typedef   MapMotClef::iterator iterator;

  public: 
  int linenumber,charnumber;
  list<string> ffincludedir;
  typedef  list<string>::iterator Iffincludedir;
  typedef  list<string>::const_iterator ICffincludedir;
  
  private:
  bool firsttime;
  int level;

  char buf[1024];
  int typetoken;
  bool echo;
  stack<char *> strdata;

  struct xxxx { 
    int l;
    istream * f;
    const string * filename; 
    istream * nf;

    xxxx() : l(0), f(0) , filename(0),nf(0)   {}   
    void  open(mylex *lexx,const char * ff) ;
    void  readin(mylex *lexx,const string & s,const string *name=0);
    void close() ;
  };
  
  friend struct mylex::xxxx;
  
  xxxx pilesource[100];
  istream & source() const {return  * pilesource[level].f;} 
  ostream & cout ;

  MapMotClef  MotClef;
  list<MapMacroDef> *listMacroDef;
  list<MapMacroParam> *listMacroParam;
  public:
  
  mylex(ostream & out);
  string token() const;
  void print(ostream &f) const; 
  int scan(int lvl=0);
  int lineno(){return linenumber;}
  char * YYText() { return buf;}
  void dump(ostream & f ) ;
  
  void erreur(const char * s) {
    cerr << " Error line number" <<linenumber << ": " << s << endl;
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
      return  pilesource[level].filename ? pilesource[level].filename->c_str() : " -- in macro -- ";
    return "-- unkown --";}
    
  void input(const char *  filename) ;
  void input(const string &str,const string *name=0);
  bool close() ;

  char * newcopy(const char * s)
   {
    char *r(new char  [strlen(s)+1]);
    strcpy(r, s);
    strdata.push(r);
    return r;
  }
 ostream & ShowStack(ostream & f); 
~mylex();
private: 
  int basescan();  
  int EatCommentAndSpace(string *data=0);
  int scan1();
  bool SetMacro(int &ret);
  bool CallMacro(int &ret);
  char * match(int i);
  void ErrorScan(const char * s) {
      cerr  << "\n" ;
      ShowStack(cerr);
    throw(ErrorCompile(s,lineno(),YYText() ) );}    
  
  
} ;

 mylex * Newlex(  ostream & out);
 void Destroylex(mylex * m);

extern mylex *zzzfff;
#endif
