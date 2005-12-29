#include  <complex>
#include  <string>
#include  <iostream>
#include  "error.hpp"
#include  <ctype.h>
#include  <stdlib.h>
#include  <map>
#include "AFunction.hpp"
//class pfes;
#include  <iomanip>                                                                                                        
#include "lg.tab.hpp"
#include "lex.hpp"


extern YYSTYPE *plglval;

//  New version of macro expantion more classical
//  and more  simple
// FH Jan. 2005

static const bool debugmacro = false;

/*inline char * newcopy(const char * s)
{
  char *r(new char  [strlen(s)+1]);
  strcpy(r, s);return r;
}
*/
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
    ffassert(r);
    return true;}}

void  mylex::Add(Key k,int r,aType t){ 
  iterator ii= MotClef.find(k);
  ffassert(ii==MotClef.end());
  MotClef.insert(make_pair<const Key,Value>(k,make_pair(r,t))); }

void mylex::dump(ostream & f ) 
{
  const_iterator i=MotClef.begin(),e=MotClef.end();
  for (;i != e;i++)
    {
      f << "      " << i->first << " " << i->second.first << " " <<  i->second.second->name() << endl;
    }
}   

int mylex::EatCommentAndSpace(string *data)
{
 // if data is !0 then add reading char in data
 // return the last read char c
 //  --------------------
  int c,caux;
  int incomment =0;
  do {
    incomment = 0;
    c = source().peek();
    
    // eat spaces 
    while (isspace(c) || c == 10 || c == 13 )
     {
      c = source().get();
      if(c==10 || c==13) c='\n';
      if (echo) cout << (char) c;
      if(c=='\n') { linenumber++; if (echo) cout << setw(5) <<linenumber << " : " ;};             
       if(data) *data+=char(c);
       c=source().peek();
      } 
     
    // eat comment  
    if(c=='/') 
      { c = source().get();
        caux=source().peek(); 
        if(caux =='/') incomment = 1;
        else if (caux == '*' ) incomment = 2;
        if(!incomment) source().putback(c);
      }
      
      
    if(incomment==1) 
      {
        if (echo) cout << "//" ;source().get();
        if(data) *data+="//";

        do {c=source().get();            
        if(c==10 || c==13) c='\n';
        if (echo) cout << (char) c;
        if(data) *data+=char(c);
        if(c=='\n') { linenumber++; if (echo) cout << setw(5) <<linenumber << " : " ;};             
        }
      while( (c!= '\n') && (c!= 10)  && (c!= 13)  &&  ( c != EOF) );
      }
    else if(incomment==2) 
      {
        if (echo) cout << "/*" ;
        if(data) *data+="/*";

        source().get();
        do {    
          c=source().get(); 
          if(c==10 || c==13) c='\n';
          if (echo) cout << (char) c ;
          if(data) *data+=char(c);
          if(c=='\n') { linenumber++; if (echo) cout << setw(5) <<linenumber << " : " ;};             
          caux = source().peek();
        } while(c != EOF && !(c=='*' && caux=='/') ) ;
        
        if(c != EOF)  
          {     c = source().get();
          if (echo) cout << (char) c ;
          if(data) *data+=char(c);          
          }
        else erreur( " Unterminated comment");
      }
  } while (incomment);
  return c;
}
int mylex::basescan()
{
  extern long mpirank;
   
  int c;
  buf[0]=0;
  buf[1]=0;
  buf[2]=0;
  buf[3]=0; //
 debut:
  TheCurrentLine=linenumber;
  // modif FH 
  if (firsttime) 
    firsttime=false,cout << setw(5) <<linenumber << " : " ; 

  EatCommentAndSpace(); 
  c =source().get(); // the current char 
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
      { buf[i++]=source().get(),buf[i]=0,plglval->dnum = atof(buf),ret=CNUM;}
    else 
      if (real)
        plglval->dnum = atof(buf),ret=DNUM;
      else
        plglval->lnum = atoi(buf),ret=LNUM;
    
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
      
    }
  else  if (c=='"')
    {
      int i;
      char cc;
      for (     i = 0,cc=source().peek(); 
                i < 256 &&  (isprint(cc) && cc !='\n'  && cc !='"');
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
              case 10:
              case 13:
		cc='\n';
		linenumber++;
		//if (echo) cout << setw(5) <<linenumber << " : " ;		
              default : buf[i]=cc  ;break;              
              }
            }
          else 
            buf[i] = source().get();
        }
      if (i == 256) erreur ("String too long");
      buf[i] = 0;
      if(source().get() != '"') erreur("End of String could not be found");
      plglval->str = newcopy(buf);
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
      case '#': 
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
//      case '_': 
      case '?':
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
      strcpy(plglval->oper,buf);
      if(lexdebug)  cout << "Special '" <<  plglval->oper << "' " << ret << " ";
    }
  typetoken=ret; 
  return ret;
}

int mylex::scan1()
{

  extern long mpirank;
 // bool echo = mpirank == 0; 

  int ret= basescan();

  if ( ret == ID)
    while ((SetMacro(ret)))0;

  if ( ret == ID)
    while ((CallMacro(ret)))0;
  return ret;
}    

int mylex::scan(int lvl)
{
 
  int ret= scan1(); 
  
      
  if ( ret == ID) {
    if (! InMotClef(plglval->type,ret))  {
      if (  FindType(buf) == 1)  
         ret =  FESPACE;
     plglval->str = newcopy(buf);
      }}
  
  if ( ret =='{') { //cout << " listMacroDef->push_back"<< endl; 
      listMacroDef->push_back( MapMacroDef() );}
  else if (ret == '}') {//cout << " listMacroDef->pop_back"<< endl;
       listMacroDef->pop_back( );}
  
  if (! lexdebug && echo && lvl==0 ) print(cout);
  
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
   int k=0;
    if (typetoken == STRING) 
      {
	
        f <<'"';
        while (buf[++i]) { k++;
	if (buf[i]=='\n') k++, f << "\\n"; 
	else if (buf[i]=='"') k++,f << "\\\""; 
	else f << buf[i] ;
	if ( k%50==49) f << "\n  ... : ";  
	}
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
	      if (rr != ID) {
		cout <<" Erreur waiting of an ID: " << buf << " " << rr <<  endl;
		erreur("in macro def: waiting for a ID"); }
	      macroparm.push_back(buf);
	      rr =  basescan(); 
	      if (rr  == ')') break;          
	      if ( rr != ',') erreur("in macro def: waiting  , or  ) "); 
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
  // New Version with direct macro expansion
  // FH  jan 2005 
  // -----------------------------------------
  for (list<MapMacroDef>::const_iterator i=listMacroDef->begin(); i != listMacroDef->end(); i++)
    {
      MapMacroDef::const_iterator j= i->find(buf);
      
      if (j != i->end()) {
	// int inmacroold=withmacropara;
	const deque<string>  & macroparm= j->second;
	int nbparam= macroparm.size()-1;
	MapMacroParam  lp;
	if (nbparam > 0 ) { 
	  match('(');
	  for (int k=0;k<nbparam;k++)
	    { 
	      string p;
	      int kend= ( k+1 == nbparam) ? ')' : ',';
	      int lvl=0;
	      while (1) {
		int rr = basescan();
		if(lvl && rr==')') lvl--; //   if ( then  we eat next ) 
		else if (rr=='(') lvl++ ;  //  eat next 		
		else if (lvl<=0) {
		  if (rr==kend ) break;
		  else if  (rr==')' || rr==',')  {// Correction FH 2/06/2004
		  cerr << "Error in macro expantion "<< j->first 
		       << ", we wait for "<< char(kend) << " and we get  " << char(rr)<< endl;
		  cerr << " number of macro parameter in definition is " << nbparam << endl;
		  ErrorScan(" Wrong number of parameter in  macro call");
		}}
		
		if (rr==ENDOFFILE) ErrorScan(" ENDOFFILE in macro usage");
		p += token(); // Correction FH 2/06/2004 of string parameter
	      }
	      if(debugmacro)
		cout << "macro arg "<< k << " :" << macroparm[k] << " -> " <<  p << endl;
	      lp.insert(make_pair<string,string>(macroparm[k],p));
	      //lp[macroparm[k]] = p; 
	    }
	}
	if(debugmacro)
	  cout <<   " input in : -> " << macroparm[nbparam]  << " " << nbparam << endl;
	input(macroparm[nbparam]);
	// ici il faut faire la substitution  de parameter 
	// -----------------------------------------------
	string expandtxt;
	bool echosave=echo;
	while(1) {
	  int c= EatCommentAndSpace(&expandtxt); // eat comment to save it;
	  if (c == EOF) break;
	  ret = basescan();
	  if(ret== ENDOFFILE) break; 
	  if (nbparam && ret == ID) 
	    {  
	      MapMacroParam::const_iterator j=lp.find(buf);
	      if ( j !=  lp.end()) 
		expandtxt+=j->second;
	      else 
		expandtxt+=buf; 
	    }
	  else if (ret!='#')  //  macro concatenation operator 
	    expandtxt+=token();
	}
	echo=echosave;
	input(expandtxt);
	ret =  scan1(); // Correction FH 6/06/2004 of string parameter
	return true;        
      }
    }
  return false;
}

void  mylex::xxxx::open(mylex *lex,const char * ff) 
{
  //filename=ff;
  l=0;
  nf=f=0;
  if(lex->ffincludedir.empty()) // if no liste 
    nf=f= new ifstream(ff); 
  if (!f || !*f)
   {
   if ( f)  delete f; 
   for (ICffincludedir k=lex->ffincludedir.begin();
        k!=lex->ffincludedir.end();
        ++k)
   {
    string dif_ff(*k+ff);
    nf=f= new ifstream(dif_ff.c_str()); 
    if ( f)  {
      if ( f->good()) {  
        filename = new string(dif_ff);
        break;
      }
      delete f;
     }     
   } 
   } 
   else
      filename=new  string(ff);
  if (!f || !*f) {
    lex->cout << " Error openning file " <<ff<< endl;
    lgerror("lex: Error input openning file ");};
}
void  mylex::xxxx::readin(mylex *lex,const string & s)//,int nbparam,int bstackparam) 
{
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
  if (filename) delete filename;
  
}
void mylex::input(const char *  filename) 
{
  ffassert(level<99 && level >= -1);
  if (level>=0) 
    pilesource[level].l =linenumber;
  
  pilesource[level+1].open(this,filename);
  pilesource[level+1].l =0;
  // cout << "\n ++include " << filename << ";" << level+1 << endl;
  linenumber = 1;     
  level++;      
}

void mylex::input(const string & str) 
{
  
  ffassert(level<99 && level >= -1);
  if (level>=0) 
    { pilesource[level].l =linenumber;
    }
  
  pilesource[level+1].readin(this,str);
  linenumber = 0;     
  level++;
  
}
  
bool mylex::close() { 
  if(debugmacro )
    cout << "\n close " << level ;
  ffassert(level >=0 && level < 100);
  // cout << "\n-- close " << level << endl;
  pilesource[level].close(); 
  // cout << "\n ++   " << level << endl;
  if (--level<0)
    return false;
  linenumber = pilesource[level].l;
  return true;
}
  
 mylex::~mylex()
{
  delete listMacroDef;
  while( ! strdata.empty()) 
   {  //  commente july 2005 FH ???. a test plus finement. 
      // delete [] strdata.top(); // bug  ????? FH  25032005 
      strdata.pop();
   }
}
