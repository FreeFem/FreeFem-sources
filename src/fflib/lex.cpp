#include  <complex>
#include  <string>
#include  <iostream>
#include  <iomanip>
#include  "error.hpp"
#include  <ctype.h>
#include  <stdlib.h>
#include  <map>
#include "AFunction.hpp"
//class pfes;
#include "lg.tab.hpp"
#include "lex.hpp"


static const bool debugmacro = false;

inline char * newcopy(const char * s)
{
  char *r(new char  [strlen(s)+1]);
  strcpy(r, s);return r;
}

void  mylex::Add(Key k,int i)   
{
  Check(!i,k,"mot clef");
  Add(k,i,0); 
}
void  mylex::Add(Key k,aType t)   
{
  Check(!t,k,"type");
  Add(k,TYPE,t);
}
void  mylex::AddF(Key k,aType t)  
{
  Check(!t,k,"FUNCTION");
  Add(k,FUNCTION,t); 
}
void  mylex::AddSpace(Key k,aType t)  
{
  Check(!t,k,"FESPACE");
  Add(k,FESPACE,t); 
}

bool mylex::InMotClef  (aType & t, int & r) const {
  const_iterator i= MotClef.find(buf);
  if (i== MotClef.end()) {
    t=0;r=ID;
    return false;}
  else {
    r=i->second.first;
    t=i->second.second;
    throwassert(r);
    return true;}}

void  mylex::Add(Key k,int r,aType t){ 
  iterator ii= MotClef.find(k);
  throwassert(ii==MotClef.end());
  MotClef.insert(make_pair<const Key,Value>(k,make_pair(r,t))); }

void mylex::dump(ostream & f ) 
{
  const_iterator i=MotClef.begin(),e=MotClef.end();
  for (;i != e;i++)
    {
      f << "      " << i->first << " " << i->second.first << " " <<  i->second.second->name() << endl;
    }
}   

int mylex::basescan()
{
  extern long mpirank;
  bool echo = mpirank == 0; 
  int c,caux;
  int incomment =0;
  buf[0]=0;
  buf[1]=0;
  buf[2]=0;
  buf[3]=0; //
 debut:
  TheCurrentLine=linenumber;
  do {
    incomment = 0;
    c = source().get();
    
    while (isspace(c) || c == 10 || c == 13 ){
      if(c==10 || c==13) c='\n';
      if (echo) cout << (char) c;
      if(c=='\n') { linenumber++; if (echo) cout << setw(5) <<linenumber << " : " ;};             
      c=source().get();} 
    if(c=='/') 
      { 
        caux=source().peek(); 
        if(caux =='/') incomment = 1;
        else if (caux == '*' ) incomment = 2;
      }
    if(incomment==1) 
      {
        if (echo) cout << "//" ;source().get();
        do {c=source().get();            
        if(c==10 || c==13) c='\n';
        if (echo) cout << (char) c;
        if(c=='\n') { linenumber++; if (echo) cout << setw(5) <<linenumber << " : " ;};             
        }
      while( (c!= '\n') && (c!= 10)  && (c!= 13)  &&  ( c != EOF) );
      }
    else if(incomment==2) 
      {
        if (echo) cout << "/*" ;
        source().get();
        do {    
          c=source().get(); 
          if(c==10 || c==13) c='\n';
          if (echo) cout << (char) c ;
          if(c=='\n') { linenumber++; if (echo) cout << setw(5) <<linenumber << " : " ;};             
          caux = source().peek();
        } while(c != EOF && !(c=='*' && caux=='/') ) ;
        
        if(c != EOF)  
          {     c = source().get();
          if (echo) cout << (char) c ;
          }
        else erreur( " Unterminated comment");
      }
  } while (incomment);
  
  char nc = source().peek(); // next char
  buf[0]=c;
  buf[1]=nc;
  buf[2]=0;
  int ret = c;
  if (c == EOF) 
    { 
      //if (echo) cout << "ENDOFFILE "<< endl;
      if (close() )  goto debut; 
      buf[0]=0;
      return ENDOFFILE;}
  else if (isdigit(c) || c=='.' && isdigit(nc)) {
    //  a number
    int i=1;
    buf[0]=c;
    bool real= (c=='.');
    
    
    while ( isdigit(nc) && i< 50 ) buf[i++]=source().get(),nc=source().peek();
    if (!real && (nc == '.')) real=true,buf[i++]=source().get(),nc=source().peek();
    while ( isdigit(nc) && i< 50 ) buf[i++]=source().get(),nc=source().peek();
    if (nc =='E' || nc =='e' || nc =='D' || nc =='d') {real=true;
    buf[i++]=source().get(),nc=source().peek();
    if (nc =='+' || nc =='-' || isdigit(nc))  buf[i++]=source().get(),nc=source().peek();
    while ( isdigit(nc) && i< 50 ) buf[i++]=source().get(),nc=source().peek();
    }
    if (i>= 50)  erreur("Number to long");
    
    buf[i]=0;
    if (nc=='i' ) 
      { buf[i++]=source().get(),buf[i]=0,lglval.dnum = atof(buf),ret=CNUM;}
    else 
      if (real)
        lglval.dnum = atof(buf),ret=DNUM;
      else
        lglval.lnum = atoi(buf),ret=LNUM;
    
    if(lexdebug) cout << " NUM : " << buf  << " "<<endl;
  }
  else if (isalpha(c) || c=='\\')
    {
      ret =  ID;
      int i; 
      for (i = 1; i < 256 && isalnum(source().peek()); i++) 
        buf[i]=source().get();
      if (i == 256) 
        erreur ("Identifyier too long");
      buf[i] = 0;
      /*      if (! InMotClef(lglval.type,ret))  {
              if (  FindType(buf) == 1)  
              ret =  FESPACE;
              lglval.str = newcopy(buf); }*/
      
    }
  else  if (c=='"')
    {
      int i;
      char cc;
      for (     i = 0,cc=source().peek(); 
                i < 256 &&  (isprint(cc) && cc !='\n' && cc !='"');
                cc=source().peek(),++i
                ) 
        {       
          if ( cc == '\\')  //   escape 
            {
              cc= source().get();
              cc= source().get();      
              switch (cc) {
              case 'n': buf[i]='\n';break;
              case 'f': buf[i]='\f';break;
              default : buf[i]=cc  ;break;              
              }
            }
          else 
            buf[i] = source().get();
        }
      if (i == 256) erreur ("String too long");
      buf[i] = 0;
      if(source().get() != '"') erreur("End of String could not be found");
      lglval.str = newcopy(buf);
      if(lexdebug)  cout << "STRING=" << '"' << buf << '"' << endl;
      
      ret= STRING;
    }
  else 
    {
      ret = c;
      switch(c) {
        
      case '{': 
      case '%': 
      case '}': 
      case '(': 
      case ')': 
      case '[': 
      case ']': 
      case ',': 
      case ';': 
        break;
      case '*': 
        if (nc == '=') ret=MULEQ;       
        break;
      case '/': 
        if (nc == '=') ret=DIVEQ;
        break;  
      case '^': 
      case '~': 
      case '\'': 
      case '_': 
        break;
      case '.': 
        if (nc == '*') ret = DOTSTAR;
        else if (nc == '/') ret = DOTSLASH;
        break;
      case ':': 
        if (nc == '=') ret= SET;
        break;
      case '&':
        if (nc == '&') ret= AND;
        break;
      case '|':
        if (nc == '|') ret= OR;
        break;
      case '!': 
        if (nc == '=') ret=NE;
        break;          
      case '<':
        if (nc == '=') ret=LE;
        if (nc == '<') ret=LTLT;
        if (nc == '>') ret= NE;
        break;
      case '>':
        if (nc == '=') ret=GE;  
        if (nc == '>') ret=GTGT;
        break;  
      case '=': 
        if (nc == '=') 
          ret=EQ;
        break;
      case '-': 
        if (nc == '>') ret=ARROW;       
        if (nc == '-') ret=MOINSMOINS;  
        if (nc == '=') ret=MOINSEQ;     
        break;
      case '+': 
        if (nc == '+') ret=PLUSPLUS;    
        if (nc == '=') ret=PLUSEQ;      
        break;
      default: 
        cerr << "'" << (char) c << (char) nc << "' <=> " << (int) c << " is " ;
        erreur (" Unexpected character");
      }
      
      if (ret!=c) source().get() ;
      else buf[1] = 0;
      strcpy(lglval.oper,buf);
      if(lexdebug)  cout << "Special '" <<  lglval.oper << "' " << ret << " ";
    }
  typetoken=ret; 
  return ret;
}
int mylex::scan()
{

  extern long mpirank;
  bool echo = mpirank == 0; 

  int ret= basescan();
  if ( ret == ID)
    while ((SetMacro(ret)))0;

  if (withmacropara)
    if (ret==ID)
       (ExpandParam(ret));
  
  if ( ret == ID)
    while ((CallMacro(ret)))0;
    
      
  if ( ret == ID) {
    if (! InMotClef(lglval.type,ret))  {
      if (  FindType(buf) == 1)  
         ret =  FESPACE;
     lglval.str = newcopy(buf); }}
  
  if ( ret =='{') { //cout << " listMacroDef->push_back"<< endl; 
      listMacroDef->push_back( MapMacroDef() );}
  else if (ret == '}') {//cout << " listMacroDef->pop_back"<< endl;
       listMacroDef->pop_back( );}
  
  if (! lexdebug && echo ) print(cout);
  
  return ret;
}
  string mylex::token() const
  {
   int i=-1;
   string f;
    if (typetoken == STRING) 
      {
        f += '"';
        while (buf[++i]) 
          if (buf[i]=='\n') f +=  "\\n"; 
          else if (buf[i]=='"') f += "\\\""; 
          else f +=  buf[i] ;
        f+= '"';
      }
    else 
      while (buf[++i]) 
        if (buf[i]=='\n') f += "\\n"; 
        else f += buf[i] ;
    return f;
  }
  
  void mylex::print(ostream &f) const
  {
   int i=-1;
    if (typetoken == STRING) 
      {
        f <<'"';
        while (buf[++i]) 
          if (buf[i]=='\n') f << "\\n"; 
          else if (buf[i]=='"') f << "\\\""; 
          else f << buf[i] ;
        f << '"';
      }
    else 
      while (buf[++i]) 
        if (buf[i]=='\n') f << "\\n"; 
        else f << buf[i] ;
  
  } 

char * mylex::match(int i) 
{ if ( i != basescan() ) {
  cerr << " we waiting for a '" << char(i) <<"'" << endl;
  ErrorScan(" err macro ");};
 return buf; }

bool mylex::ExpandParam(int &ret)
{
  if (listMacroParam )
  {
    string arg=buf;
    for(list<MapMacroParam>::const_iterator i=listMacroParam->begin(); i != listMacroParam->end(); i++)
      {
        MapMacroParam::const_iterator j=i->find(arg);
        if ( j != i->end() )
          {
            if(debugmacro)
             cout << "macro arg : " <<  arg << " -> " << j->second << endl;
            input(j->second, 0);
            ret = basescan(); 
            arg=buf;
            //return true;
          }     
      }
  }
      return false;
      }
  
  bool mylex::SetMacro(int &ret)
    { 
      bool rt=false;
      if (strcmp(buf,"macro")==0)
        {
          char *macroname=newcopy(match(ID));
          int nbparam =0;
          deque<string> macroparm;
          int rr=basescan(); 
          string def;                          
          if (rr!='(') 
           def += buf;
           else
            { // a '(' after macro name
            rr =  basescan();
            
            if (rr != ')' ) 
              do { 
                if (nbparam++>= 100) 
                  { cerr << "in macro " <<macroname << endl;
                  ErrorScan(" Too much (more than 100) parameters"); }  
                
                //cout << " \t\t " << ID << " " << rr << " " << buf << endl;
                if (rr != ID)  
                  { cout <<" Erreur waiting of an ID: " << buf << " " << rr <<  endl;
                    erreur("in macro def: waiting for a ID"); }
                macroparm.push_back(buf);
                rr =  basescan(); 
                if (rr  == ')') break;          
                if ( rr != ',') erreur("in macro def: waiting for , or ) "); 
                rr =  basescan();
                
              } while(1); 
          }
                
          do {
            int i = source().get();
            if (i == ENDOFFILE) {  cerr << "in macro " <<macroname <<  endl;
            ErrorScan(" ENDOFFILE in macro definition. remark:a macro end with // ");}
            int ii = source().peek();
            if (i == '/' && ii == '/') { source().putback('/'); break;} 
            def +=char(i);        
          } while(1);
          macroparm.push_back(def);
          cout << "macro " << macroname  ;
          for (int i=0;i<macroparm.size()-1;i++)
            cout << ( (i == 0) ? '(' : ',') << macroparm[i];
          if (nbparam)
            cout << " )  " ;
          for (int i=0;i<def.size();i++)
            if (def[i] == 10 || (def[i] == 13 ) ) 
             {
              def[i]='\n';         
             linenumber++; cout << '\n' << setw(5) <<linenumber << " : " ;
             }            

            else 
              cout << def[i]   ;
         // cout << macroparm.back() ;
          MapMacroDef & MacroDef =listMacroDef->back();
          MapMacroDef::const_iterator i=MacroDef.find(macroname);
          if ( i == MacroDef.end() )
            MacroDef[macroname]=macroparm;
          else { 
            cerr << "ERREUR in macro name:" <<macroname << endl;
            ErrorScan(" The macro already exist, sorry");}
          rt=true;
          ret= basescan();
        }
      return rt;
    }
  bool mylex::CallMacro(int &ret)
    {
      
      for (list<MapMacroDef>::const_iterator i=listMacroDef->begin(); i != listMacroDef->end(); i++)
        {
          MapMacroDef::const_iterator j= i->find(buf);
          
          if (j != i->end()) {
            bool inmacroold=withmacropara;
            const deque<string>  & macroparm= j->second;
            int nbparam= macroparm.size()-1;
            if (nbparam > 0 ) { 
              match('(');
              if (!listMacroParam) listMacroParam= new  list<MapMacroParam>;
              listMacroParam->push_front(MapMacroParam());
              MapMacroParam & lp=listMacroParam->front();
              for (int k=0;k<nbparam;k++)
                { 
                  string p;
                  int kend= ( k+1 == nbparam) ? ')' : ',';
                  while (1) {
                    int rr = basescan();
                    if (rr==kend) break;
                    if (rr==ENDOFFILE) ErrorScan(" ENDOFFILE in macro usage");
                    p += buf;
                  }
                  if(debugmacro)
                  cout << "macro arg "<< k << " :" << macroparm[k] << " -> " <<  p << endl;
                  lp.insert(make_pair<string,string>(macroparm[k],p));
                  //lp[macroparm[k]] = p; 
                }
            }
            if(debugmacro)
            cout <<   " input in : -> " << macroparm[nbparam]  << " " << nbparam << endl;
            input(macroparm[nbparam], nbparam);
            ret =  basescan();
            return true;        
          }
        }
      return false;
    }
  
  void  mylex::xxxx::open(mylex *lex,const char * ff) 
    {
      filename=ff;
      l=0;
      nf=f= new ifstream(ff);     
      if (!f || !*f) {
        lex->cout << " Error openning file " <<filename<< endl;
        lgerror("lex: Error input openning file ");};
    }
  void  mylex::xxxx::readin(mylex *lex,const string & s,int nbparam) 
    {
      cas= nbparam;
      filename=0;
      l=0;
      nf=f= new istringstream(s.c_str()); 
      
      if (!f || !*f) {
        lex->cout << " Error readin string  :" <<s<< endl;
        lgerror("lex: Error readin macro ");};
    }
    
  void mylex::xxxx::close() 
    { 
      if( nf)  delete nf;
      if (filename) delete [] filename;
      
    }
  void mylex::input(const char *  filename) 
    {
      assert(level<99 && level >= -1);
      if (level>=0) 
        pilesource[level].l =linenumber;
        
      pilesource[level+1].open(this,filename);
      pilesource[level+1].cas =0;
      pilesource[level+1].l =0;
     // cout << "\n ++include " << filename << ";" << level+1 << endl;
      linenumber = 0;     
      level++;
      
    }
  void mylex::input(const string & str,int nbparam) 
    {
      withmacropara= nbparam>0; ; 
      assert(level<99 && level >= -1);
      if (level>=0) 
        { pilesource[level].l =linenumber;
        }
      
      pilesource[level+1].readin(this,str,nbparam);
      linenumber = 0;     
      level++;
      
    }
  
  bool mylex::close() { 
      if(debugmacro)
       cout << "\n close " << level << " " << pilesource[level].cas << " " << (listMacroParam ==0 ? 0 :listMacroParam->size()) << endl;
      assert(level >=0 && level < 100);
      if (pilesource[level].cas>0) //  macro with parameter 
        {
          assert(listMacroParam);
          listMacroParam->pop_front();
           if(debugmacro) cout << " " << listMacroParam->size();
         }
       if(debugmacro) cout << endl;
     // cout << "\n-- close " << level << endl;
      pilesource[level].close(); 
     // cout << "\n ++   " << level << endl;
      if (--level<0)
        return false;
      linenumber = pilesource[level].l;
       withmacropara = listMacroParam && listMacroParam->size();
       return true;}
  
