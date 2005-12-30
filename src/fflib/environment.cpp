#include "environment.hpp"
#include "iostream"
#ifdef HAVE_GETENV
#include <stdlib.h>
#endif
using namespace std;

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
   EnvironmentData::iterator ekey=environment.find(key);
   if( ekey != environment.end()) 
    {
     OneEnvironmentData * pl= &ekey->second;
     OneEnvironmentData::iterator i=find(pl->begin(),pl->end(),item);
       return i != pl->end();
     }
    
   return false;
 }
 
bool EnvironmentInsert(string key,string item,string before)
{
   bool ret=true;
   OneEnvironmentData  & l = environment[key];
   
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

int GetEnvironment(const char * key,string items)
{
  if(verbosity>=100)  cout << key << " -> " << items << endl;
  int d=0, k=0,i;
  items+=";;";
  for ( i=0;i<items.size();i++)
   if(  items[i]==';')
    { 
      string item=TransDir(items.substr(d,i-d));
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

void GetEnvironment()
 {
 char  * ff_verbosity=0,* ff_loadpath=0,* ff_incpath=0;
#ifdef HAVE_GETENV
  ff_verbosity = getenv("FF_VERBOSITY");
  ff_loadpath = getenv("FF_LOADPATH");
  ff_incpath = getenv("FF_INCLUDEPATH");
#endif

#ifdef PURE_WIN32 
 
  const int LEN = 4096;
  char envv[LEN];
  char envl[LEN];
  char envi[LEN];

  if (GetEnvironmentVariable("FF_VERBOSITY", envv, LEN) > 0) 
   ff_verbosity=envv;
   
  if (GetEnvironmentVariable("FF_LOADPATH", envl, LEN) > 0) 
   ff_verbosity=envl;
   
  if (GetEnvironmentVariable("FF_INCLUDEPATH", envi, LEN) > 0) 
   ff_verbosity=envi;

#endif 

 if ( ff_verbosity ) { 
        verbosity = atoi(ff_verbosity);
        if(verbosity>2) cout << " --  verbosity is set to " << verbosity << endl;
  }
 if(ff_loadpath)
    GetEnvironment("load",ff_loadpath);
 if(ff_incpath)
    GetEnvironment("include",ff_incpath);
 if( verbosity >2) 
   {
    EnvironmentData::iterator load=environment.find("load");
    EnvironmentData::iterator inc=environment.find("include");    
    if(  load != environment.end()) {
      show("\nload path : ",load->second, "\n \t ");
      cout <<"(.)"<<endl;
    }
    if(  inc != environment.end()) {
      show("\ninclude path : ",inc->second, "\n \t ");
      cout <<"(.)"<<endl;}
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
