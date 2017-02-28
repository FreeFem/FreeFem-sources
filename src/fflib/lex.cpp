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

// <<mylex_Add_Key_aType>>
void  mylex::Add(Key k,aType t)   
{
  Check(!t,k,"type");

  // TYPE defined at [[file:../lglib/lg.ypp::token type TYPE]]
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

// <<mylex_InMotClef>> looks in MotClef for the token contained in buf (buf is set by by [[mylex_basescan]]) and returns
// its type as t and the grammar token id as r (eg [[file:~/ff/src/lglib/lg.ypp::TYPE]]).

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

// alh - 5/4/15 - <<mylex::InMotClef_string>> [[file:lex.hpp::InMotClef_string]] Same as [[mylex_InMotClef]], but checks
// another buffer than buf
bool mylex::InMotClef(const char *b,aType &t,int &r)const{
  const_iterator i=MotClef.find(b);
  if(i==MotClef.end()){t=0;r=ID;return false;}
  else {
    r=i->second.first;
    t=i->second.second;
    ffassert(r);
    return true;
  }
}

// <<mylex_Add_Key_int_aType>>
void  mylex::Add(Key k,int r,aType t){ 
  iterator ii= MotClef.find(k);
  if(ii!=MotClef.end())
  {
    cout << "erreur add motclef " <<  k<< " " << r << " " << *t << endl; 
    ffassert(ii==MotClef.end());
  }
  //ffassert(ii==MotClef.end());
  MotClef.insert(pair<const Key,Value>(k,Value(r,t))); }

// <<dump>>
void mylex::dump(ostream & f ) 
{
  const_iterator i=MotClef.begin(),e=MotClef.end();
  for (;i != e;i++)
    {
      f << "      " << i->first << " " << i->second.first << " " <<  i->second.second->name() << endl;
    }
}   

inline bool isNLCR(istream & f,int c)
{
    // eat  CR NL or NL CR paire  
    int cc= f.peek();
    bool ret=(c == 10 || c == 13) ;
    if(ret && ( cc != c) && (cc == 10 || cc == 13) )
     f.get();
    return ret;
}

int mylex::EatCommentAndSpace(string *data)
{
 // if data is !0 then add reading char in data
 // return the last read char c
 //  --------------------
  int c,caux,sep=0;
    const int space=(int) ' ';
  int incomment =0;
  do {
    incomment = 0;
    c = source().peek();
    
    // eat spaces 
    while (isspace(c) || c == 10 || c == 13 )
      {sep=space;
      c = source().get();
      if(isNLCR(source(),c)) c='\n';
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
      { sep=space;
        if (echo) cout << "//" ;source().get();
        if(data) *data+="//";

        do {c=source().get();            
        if(isNLCR(source(),c)) c='\n';
        if (echo) cout << (char) c;
        if(data) *data+=char(c);
        if(c=='\n') { linenumber++; if (echo) cout << setw(5) <<linenumber << " : " ;};             
        }
      while( (c!= '\n') && (c!= 10)  && (c!= 13)  &&  ( c != EOF) );
      }
    else if(incomment==2) 
      { sep=space;
        if (echo) cout << "/*" ;
        if(data) *data+="/*";

        source().get();
        do {    
          c=source().get(); 
          if(isNLCR(source(),c)) c='\n';
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
    return (c==EOF) ? c : sep;
}

// <<mylex_basescan>>
int mylex::basescan()
{
  //  extern long mpirank;
   
  int c;
  buf[0]=0;
  buf[1]=0;
  buf[2]=0;
  buf[3]=0; //
 debut:
  TheCurrentLine=linenumber;
  // modif FH 
  if (firsttime) 
  {
      firsttime=false;
    if(echo) cout << setw(5) <<linenumber << " : " ; 
  }
  EatCommentAndSpace(); // [[mylex::EatCommentAndSpace]]
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
      return ENDOFFILE;
    }

  // <<found_a_number>>
  else if (isdigit(c) || (c=='.' && isdigit(nc))) {
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
    if (i>= 50) erreur("Number too long");
    
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

  // <<found_an_identifier>>
  else if (isalpha(c) || c=='\\')
    {
      ret =  ID;
      int i; 
      for (i = 1; i < 256 && isalnum(source().peek()); i++) 
        buf[i]=source().get();
      if (i == 256) 
        erreur ("Identifier too long");
      buf[i] = 0;
    }

  // <<found_a_string>>
  else if (c=='"')
    {
      int i;
      char cc,ccc;
      for (     i = 0,cc=source().peek(); 
                i < 256 &&  ( (isprint(cc)|isspace(cc)) && cc !='\n'  && cc !='"');
                cc=source().peek(),++i
                ) 
        {       
          if ( cc == '\\')  //   escape 
            {
              cc= source().get();
              cc= source().get(); 
	      ccc= source().peek(); 
              switch (cc) {
              case 'n': buf[i]='\n';break;
              case 'r': buf[i]='\r';break;
              case 'f': buf[i]='\f';break;
              case 't': buf[i]='\t';break;
              case 'a': buf[i]='\a';break;
              case 10:
              case 13:
		cc='\n';
		if(ccc!=cc && (ccc==10 || ccc==13)) source().get();// NL CR eat ...
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

  // Anything else (eg brackets and operators)
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
      case '_': 
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
      if( (ret == DOTSTAR) || (ret==DOTSLASH))
      {
          source().get();
          nc = source().peek();
          if(nc == '=' )
          {
              buf[2]='=';// ad FH 19 april 2012 (bug in macro ) 
              buf[3]=0; 
              source().get();
              ret = (ret == DOTSTAR) ?DOTMULEQ : DOTDIVEQ;
          }
      }
      else if (ret!=c) source().get();
      else buf[1] = 0;
      strcpy(plglval->oper,buf);
      if(lexdebug)  cout << "Special '" <<  plglval->oper << "' " << ret << " ";
    }
  typetoken=ret; 
  return ret;
}

// <<mylex_scan1>>
int mylex::scan1()
{

  // extern long mpirank;
 // bool echo = mpirank == 0; 

  int ret= basescan(); // [[mylex_basescan]]
  if(debugmacro)  cout << " scan1 " << ret << " " << token() << " " << ID << endl;
    while ( ret == ID &&SetMacro(ret)); // correction mars 2014 FH
    while ( ret == ID && CallMacro(ret)) ; // correction mars 2014 FH
    
  return ret;
}    

// <<mylex_scan>>
int mylex::scan(int lvl)
{
 
  int ret= scan1(); 
  
  // ID defined at [[file:../lglib/lg.ypp::ID]] and detected at [[found_an_identifier]]
  if ( ret == ID) {
    if (! InMotClef(plglval->type,ret))  {
      int ft = FindType(buf);

      // FESPACE, FESPACE1, FESPACE3 defined at [[file:../lglib/lg.ypp::FESPACE]]
      int feid3[4]  ={ ID,FESPACE1,FESPACE,FESPACE3};

      assert ( ft >= 0 && ft <= 3)  ;
      ret =  feid3[ft];      
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

    // found at [[found_a_string]]
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
{ if ( i != basescan() ) {// basescan -> scan1 change 2/2/2007  (non pas des substitution de parametres FH)
  cerr << " we waiting for a '" << char(i) <<"'" << endl;
  ErrorScan(" err macro ");};
 return buf; }


bool mylex::SetMacro(int &ret)
{ 
  char endmacro[]="EndMacro";
  char newmacro[]="NewMacro";
    
  bool rt=false;
  int oldmacro=1;
  if (strcmp(buf,"macro")==0 || (oldmacro=strcmp(buf,newmacro))==0 )
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
		cerr <<" Erreur waiting of an ID: " << buf << " " << rr <<  endl;
		erreur("in macro def: waiting for a ID"); }
	      macroparm.push_back(buf);
	      rr =  basescan(); 
	      if (rr  == ')') break;          
	      if ( rr != ',') erreur("in macro def: waiting  , or  ) "); 
	      rr =  basescan();
	      
	    } while(1); 
	}
      int kmacro=0;
      

      do {
	int lk=0;
	string item;
	int i = source().get();
	if (i == EOF) {  cerr << "in macro " <<macroname <<  endl;
	ErrorScan(" ENDOFFILE in macro definition. remark:a macro end with // ");}
	int ii = source().peek();
	if(isspace(i) && isalpha(ii) ) {
	    def +=char(i);
	    i = source().get();
	    item = "";
	    while(isalpha(i))
	      {
		item += char(i);
		i = source().get();
	      }
	    if( item == newmacro)  kmacro++;
	    if( item == endmacro)  {
		if (kmacro==0)  
		    { source().putback(i); break;}
		kmacro--;
	    }
	    def += item;
	    item ="";
	    ii = source().peek();
	}
	
	if(oldmacro)
	  {
	    if (i == '/' && ii == '/') { source().putback('/'); break;} 
	  }

	def +=char(i);        
      } while(1);
      macroparm.push_back(def);
      if(echo) cout << "macro " << macroname  ;
      for (size_t i=0;i<macroparm.size()-1;i++)
	if(echo) cout << ( (i == 0) ? '(' : ',') << macroparm[i];
      if (nbparam)
	if(echo) cout << " )  " ;
      for (size_t i=0;i<def.size();i++)
	if (def[i] == 10 || (def[i] == 13 ) ) 
	  {
	    def[i]='\n';         
	    linenumber++; if(echo) cout << '\n' << setw(5) <<linenumber << " : " ;
	  }                  
	else 
	  if(echo) cout << def[i]   ;
         // cout << macroparm.back() ;
      MapMacroDef & MacroDef =listMacroDef->back();
      MapMacroDef::const_iterator i=MacroDef.find(macroname);
      if ( i == MacroDef.end() )
	MacroDef[macroname]=macroparm;
      else { 
	cerr << "ERREUR in macro name:" <<macroname << endl;
	ErrorScan(" The macro already exists, sorry");}
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
//  add Stringification,FILE, LINE  march 2014 FH..
  if(strcmp(buf,"Stringification")==0)
  {
      if(debugmacro) cout <<"call Stringification : " << buf << endl;

      if(echo) cout << buf << "(" ;
      string p,cmm;
      int lvll=0;
           match('(');
      
      int kend=')';
      while (1) {
          cmm="";
          int sep =  EatCommentAndSpace(&cmm);
          p += cmm;
          int rr = scan1();// do macro expantion
          if(lvll && rr==')') lvll--; //   if ( then  we eat next )
          else if (rr=='(') lvll++ ;  //  eat next
          else if (lvll<=0)
              {
                  if (rr==kend ) break;
                  else if  (rr==')' || rr==',')  {
                      cerr << "Error in Stringification  "
                      << ", we wait for "<< char(kend) << " and we get  " << char(rr)<< endl;
                      ErrorScan(" Wrong number of parameter in  Stringification call");
                  }}
          if(debugmacro)cout << " ..." << rr << " " << token()<< " " << level << endl;
          if (rr==ENDOFFILE) ErrorScan(" Stringification in macro ");
          if(echo) cout << token();
          p += token();
          //if(echo) cout <<buf;
      }
      
      plglval->str = newcopy(p.c_str());
      ret = STRING;

      return false;
  }

  // <<FILE_macro>>
  else if(strcmp(buf,"FILE")==0)
  {
      plglval->str = newcopy(filename() );
      ret = STRING;
     return false;
  }

  // <<LINE_macro>>
  else if(strcmp(buf,"LINE")==0)
  {
    plglval->lnum = linenumber;
    ret=LNUM;
    return false;
  }
  else
  for (list<MapMacroDef>::const_iterator i=listMacroDef->begin(); i != listMacroDef->end(); i++)
    {
      MapMacroDef::const_iterator j= i->find(buf);
      
      if (j != i->end()) {
	// int inmacroold=withmacropara;
	  if(debugmacro) cout <<"call macro : " << buf << endl;
	  string * macronn= newstring(" macro: ");
	  *macronn+= buf;
	  
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
		  int lvll=0;
		  while (1) {
		      int sep =  EatCommentAndSpace();
		      int rr = basescan();// basescan -> scan1 change 2/2/2007  ( not change to pass macro name as a parameter)
		      if ( (rr=='}') && (--lvll==0) ) 
			   continue; // remove first {
		      else if ( (rr=='{') && (++lvll==1) )
			  continue; // remove last }	
		      else if(lvll==0) //  count the open close () [] 
			{  
			if(lvl && rr==')') lvl--; //   if ( then  we eat next ) 
			else if(lvl && rr==']') lvl--; //   if ( then  we eat next ) 
			else if (rr=='(') lvl++ ;  //  eat next 
			else if (rr=='[') lvl++ ;  //  eat next 
			else if (lvl<=0) 
			  {
			    if (rr==kend ) break;
			    else if  (rr==')' || rr==',')  {// Correction FH 2/06/2004
				cerr << "Error in macro expansion "<< j->first
				<< ", we wait for "<< char(kend) << " and we get  " << char(rr)<< endl;
				cerr << " number of macro parameter in definition is " << nbparam << endl;
				ErrorScan(" Wrong number of parameter in  macro call");
			    }}}
		      
		      if (rr==ENDOFFILE) ErrorScan(" ENDOFFILE in macro usage");
		      if(sep==' ') p+=' ';
		      p += token(); // Correction FH 2/06/2004 of string parameter
		      
		  }
		  if(debugmacro)
		cout << "macro arg "<< k << " :" << macroparm[k] << " -> " <<  p << endl;
	      lp.insert(pair<string,string>(macroparm[k],p));
	      //lp[macroparm[k]] = p; 
	    }
	}
	if(debugmacro)
	  cout <<   " input in : -> " << macroparm[nbparam]  << " <-> " << nbparam << endl;
	input(macroparm[nbparam]);
	// ici il faut faire la substitution  de parameter 
	// -----------------------------------------------
	string expandtxt;
	bool echosave=echo;
	while(1) {
	  int c= EatCommentAndSpace(&expandtxt); // eat comment to save it;
	  if (c == EOF) break;
	  ret = basescan();
          if(debugmacro)  cout << " ret = " << ret << token() << endl;
	  if(ret== ENDOFFILE) break; 
	  if (nbparam && ret == ID) 
	    {  
	      MapMacroParam::const_iterator j=lp.find(buf);
	      if ( j !=  lp.end()) 
		expandtxt+=j->second;
	      else 
		expandtxt+=token();
	    }
	  else if (ret!='#')  //  macro concatenation operator 
	    expandtxt+=token();
	}
	echo=echosave;
	if(debugmacro) cout <<" (macro) eadin : " << expandtxt << endl;
	input(expandtxt,macronn);
	ret =  scan1(); // Correction FH 6/06/2004 of string parameter
	return true;        
      }
    }
  return false;
}

void  mylex::xxxx::open(mylex *lex,const char * ff) 
{
  l=0;
  nf=f=0;

  // Try to open the given file name without any extra path from [[file:lex.hpp::ffincludedir]]
  if(lex->ffincludedir.empty()) // if no liste 
    nf=f= new ifstream(ff,ios_base::binary); //  modif of win32

  // If it worked, set [[file:lex.hpp::filename]] to ff, otherwise try to add path elements from
  // [[file:lex.hpp::ffincludedir]]
  if (!f || !*f)
    {
   if ( f)  { delete f; nf=f=0; }

      // for every potential path element
      for (ICffincludedir k=lex->ffincludedir.begin();
	   k!=lex->ffincludedir.end();
	   ++k)
	{
	  string dif_ff(*k+ff);
	  if (verbosity>=50) lex->cout  << "  --lex open :" << dif_ff << endl;
	  nf=f= new ifstream(dif_ff.c_str(),ios_base::binary); 
	  if ( f)  {

	    // If path works, set [[file:lex.hpp::filename]] and close test stream
	    if ( f->good()) {  
	      filename = newstring(dif_ff);
	      break;
	    }
      delete f; nf=f=0;
	  }     
	} 
    } 
  else
    filename=newstring(ff);

  // Error message if no path was right
  if (!f || !*f) {
    lex->cout << " Error openning file " <<ff<< " in: " <<endl;
    for (ICffincludedir k=lex->ffincludedir.begin();
	 k!=lex->ffincludedir.end();
	 ++k)
      lex->cout  << "  -- try  :\"" << *k+ff  << "\"\n";

    lgerror("lex: Error input openning file ");
  }
}

void mylex::xxxx::readin(mylex *lex,const string & s,const string *name, int macroargs)
{
  filename=name;
  macroarg=macroargs;
  l=0;
  nf=f= new istringstream(s.c_str()); 
  
  if (!f || !*f) {
    lex->cout << " Error readin string  :" <<s<< endl;
    if(f) delete f;
    nf = f =0;
    lgerror("lex: Error readin macro ");
  }
}

void mylex::xxxx::close() 
{ 
  // ALH_BUG Why does this segfault under Windows? Probably needs valgrind to find out.
#if !defined(ENABLE_FFCS) || !defined(WIN32)
  if(nf) delete nf;
#endif
  if(filename && (macroarg==0)) delete filename; // [[file:lex.hpp::filename]]
  f=nf=0;
}

// <<mylex_input_filename>>

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

// <<mylex_input_string>>

void mylex::input(const string & str,const string * name) 
{
  
  ffassert(level<99 && level >= -1);
  if (level>=0) 
    { pilesource[level].l =linenumber;
    }
  
  pilesource[level+1].readin(this,str,name);
  linenumber = 0;     
  level++;
  
}
  
bool mylex::close() { 
  if(debugmacro )
    cout << "\n close " << level ;
  ffassert(level >=0 && level < 100);
  // cout << "\n-- close " << level << endl;
  pilesource[level].close(); // [[mylex::xxxx::close]]
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
      delete [] strdata.top(); // bug  ????? FH  25032005 
      strdata.pop();
   }
}

 mylex::mylex(ostream & out,bool eecho):
    linenumber(1),
    charnumber(0),
    ffincludedir(ffenvironment["includepath"]),
    firsttime(true),
    level(-1),
    echo(eecho && (mpirank == 0) && verbosity),
    cout(out),
    listMacroDef(new list<MapMacroDef>),
    listMacroParam(0)
 {    
    listMacroDef->push_front(MapMacroDef());
   };
   
mylex * Newlex(  ostream & out,bool eecho)
  {
    return new mylex(out,eecho);
  }
void Destroylex(mylex * m)
 {
   delete m;
 }
ostream & mylex::ShowStack(ostream & f){
    for (int i=0;i<level;++i)
	if(  pilesource[i].filename ) f << " \t"<<i<<"\t"<<" in " << *pilesource[i].filename<< " line : " <<  pilesource[i].l << endl;
    return f;     
}
