#include "environment.hpp"
#include "iostream"
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>




#ifdef HAVE_GETENV
#include <stdlib.h>
#endif
using namespace std;
bool load(string s);

const char SLACH='/';
const char BACKSLACH='\\';

string TransDir(string dir)
{
#ifdef PURE_WIN32
  char sep=BACKSLACH, nsep=SLACH;
#else
  char nsep=BACKSLACH, sep=SLACH;
#endif
  for (int i=0; i<dir.size(); ++i)
    if(dir[i]==nsep) dir[i]=sep;
  if(dir.size()>1 && dir[dir.size()-1] != sep) 
    dir += sep;
  return  dir;    
}



template<typename T> 
void  show(char * s,const T & l,const char * separateur="\n")
{
  cout << s << * separateur;
  for (typename T::const_iterator i=l.begin(); i != l.end(); i++)
    cout  << * i << * separateur;
  //cout << endl; 
}





bool  EnvironmentFind(string key,string item)
 {
   EnvironmentData::iterator ekey=ffenvironment.find(key);
   if( ekey != ffenvironment.end()) 
    {
     OneEnvironmentData * pl= &ekey->second;
     OneEnvironmentData::iterator i=find(pl->begin(),pl->end(),item);
       return i != pl->end();
     }
    
   return false;
 }


bool  EnvironmentClean(string key)
{
   EnvironmentData::iterator ekey=ffenvironment.find(key);
   if( ekey != ffenvironment.end()) 
    {
      OneEnvironmentData * pl= &ekey->second;
      pl->clear();
      return true;
     }
    
   return false;
 }
 
bool EnvironmentInsert(string key,string item,string before)
{
   bool ret=true;
   OneEnvironmentData  & l = ffenvironment[key];
   
   OneEnvironmentData::iterator i=find(l.begin(),l.end(),item);
   
   if(i!=l.end()) {ret=false; l.erase(i);} // if existe remove
   i=find(l.begin(),l.end(),before);
   if(verbosity>=100) cout << " insert " << key << " " << item << " " << before << endl;
   if(i == l.end() && before!="$")
       l.insert(l.begin(),item); // insert in front 
   else
       l.insert(i,item); // insert before i
    
  return ret;  
}

int GetEnvironment(const string & key, string items)
{
  if(verbosity>=100)  cout << key << " -> " << items << endl;
  bool path=key.find("path")!= string::npos;
  int d=0, k=0,i;
  if(path)
    items+=";;";
  for ( i=0;i<items.size();i++)
   if(  items[i]==';')
    { 

      string item =items.substr(d,i-d);
      if(path) item=TransDir(item);
      if(verbosity>=100) cout << " + " << item << endl;
      if(!EnvironmentFind(key,item)) 
	{
	  EnvironmentInsert(key,item,"$");
	  k++; 
	}
      d=i+1;          
    }
    
 return k;   
}
int  readinitfile(const string & file)
{
	string line="";
	string key;
	string value;
	bool add=false;
	ifstream f(file.c_str());
	if( ! f) {
		if(verbosity>=20) cout << "error opening init  file: " << file << endl;
		return 0;}
	
	char c,c1,bv=0;
	int linenumber = 1;
	bool inkey=false,invalue=false,cmm=false;
	while (f)
	{
		c= f.get();
		c1=f.peek();
		//cout << c ;
		//      if(bv)
		//cerr << bv ;
		//else cerr << '.';
		if(c == EOF)
		  break;
		if(c =='\n' || c=='\r' )
		  {
		    linenumber++;
		    line="";
		    cmm = false;
		  }
		else 
		  line+= c;
		if(c=='#') cmm=true;
		if(!cmm) 
		  if(invalue) // get value  key [=|+= value ]
		    { 
		      if (bv) //  store the value
			{	  
			  if (! (c == bv || ( bv==' ' && isspace(c)) ) )
			    {
			      value+= c; // add to value
			    }
			  else  // end of value
			    {
			      bv =0; invalue=0; // fin de la value 
			      if(verbosity >= 50) 
				cout <<file <<":" << key << " = " << value <<endl;
			      if(key=="verbosity")
				verbosity=atoi(value.c_str());
			      else 
				{
				  if( !add) 
				    {EnvironmentClean(key);   
				    GetEnvironment(key,value);}
				  else 
				    {
				      bool path=key.find("path")!= string::npos;
				      if(path)
					EnvironmentInsert(key,TransDir(value),"$");		      
				      else 
					EnvironmentInsert(key,value,"$");		      
				    }
				  
				}

			      key="";
			      value="";
			    }
			}
		      else // find begin of value
			{ 
			  if(c=='\'' || c == '"' ) {bv = c;value="";}
			  else if (! isspace(c)) {value=c;bv=' ';}
			  
			}
		    }
		  else if( inkey) 
		    {
		      if ( isalnum(c) ) key += c;
		      else if ( c == '=') 
			{ 
			  inkey=false;
			  invalue=true;
			  add=false;
			}
		      else if ( c == '+' && c1=='=') 
			{
			  inkey=false;
			  invalue=true;
			  add=true;
			  c1=f.get();
			}
		      else if (! isspace(c) ) 
			break;      
		    }
		  else if(isalpha(c) ) 
		    { inkey=1; key= c;}
	}
	if( inkey || invalue || bv )
	  {
	    cout << " error read init file : " << file << " " <<  linenumber << endl;
	    cout << " line : " << line << endl;	    
	    return -11;
	  }
	return 1;
}


void GetEnvironment()
 {
 char  * ff_verbosity=0,* ff_loadpath=0,* ff_incpath=0,* home=0;
  string  ffpref="freefem++.pref";
#ifdef HAVE_GETENV
  ff_verbosity = getenv("FF_VERBOSITY");
  ff_loadpath = getenv("FF_LOADPATH");
  ff_incpath = getenv("FF_INCLUDEPATH");
  home    = getenv("HOME");
#endif

#ifdef PURE_WIN32 
 
  const int LEN = 4096;
  char envv[LEN];
  char envl[LEN];
  char envi[LEN];
  char envh[LEN];
  

  if (GetEnvironmentVariable("FF_VERBOSITY", envv, LEN) > 0) 
   ff_verbosity=envv;
   
  if (GetEnvironmentVariable("FF_LOADPATH", envl, LEN) > 0) 
   ff_verbosity=envl;
   
  if (GetEnvironmentVariable("FF_INCLUDEPATH", envi, LEN) > 0) 
   ff_verbosity=envi;
  
  if (GetEnvironmentVariable("HOMEPATH", envh, LEN) > 0) 
	  home=envh;
#endif 

#ifdef PURE_WIN32 
#else
  EnvironmentInsert("init-files","/etc/"+ffpref,"$");
#endif
  if(home) 
	EnvironmentInsert("init-files",TransDir(home)+"."+ffpref,"$");
  EnvironmentInsert("init-files",ffpref,"$");
	
  { 
  OneEnvironmentData  & l = ffenvironment["init-files"];
  OneEnvironmentData::iterator i=l.begin();	
  while( i != l.end())
    {
      if(verbosity>2) cout << " try initfile : " <<*i << endl; 
      readinitfile(*i++);	   
    }
  }
	 
	 

 if ( ff_verbosity ) { 
        verbosity = atoi(ff_verbosity);
        if(verbosity>2) cout << " --  verbosity is set to " << verbosity << endl;
  }
 if(ff_loadpath)
    GetEnvironment("loadpath",ff_loadpath);
 if(ff_incpath)
    GetEnvironment("includepaths",ff_incpath);
 if( verbosity >2) 
   {
    EnvironmentData::iterator loadpath=ffenvironment.find("loadpath");
    EnvironmentData::iterator inc=ffenvironment.find("includepath");    
    if(  loadpath != ffenvironment.end()) {
      show("\nload path : ",loadpath->second, "\n \t ");
      cout <<"(.)"<<endl;
    }
    if(  inc != ffenvironment.end()) {
      show("\ninclude path : ",inc->second, "\n \t ");
      cout <<"(.)"<<endl;}
   }
 {
   EnvironmentData::iterator toload=ffenvironment.find("load");
   if(  toload != ffenvironment.end()) 
     
     for (OneEnvironmentData::iterator i=toload->second.begin(); i != toload->second.end(); ++i)
       {
	 if(verbosity) cout << " load "<< *i << endl;
	 load(*i);
       }
   
 }
 }
 
#ifdef TESTMAIN
long verbosity=50;
 EnvironmentData  environment;
int main()
{
  GetEnvironment();
  return 0; 
}

#endif
