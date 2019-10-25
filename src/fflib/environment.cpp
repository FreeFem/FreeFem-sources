/****************************************************************************/
/* This file is part of FreeFEM.                                            */
/*                                                                          */
/* FreeFEM is free software: you can redistribute it and/or modify          */
/* it under the terms of the GNU Lesser General Public License as           */
/* published by the Free Software Foundation, either version 3 of           */
/* the License, or (at your option) any later version.                      */
/*                                                                          */
/* FreeFEM is distributed in the hope that it will be useful,               */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of           */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            */
/* GNU Lesser General Public License for more details.                      */
/*                                                                          */
/* You should have received a copy of the GNU Lesser General Public License */
/* along with FreeFEM. If not, see <http://www.gnu.org/licenses/>.          */
/****************************************************************************/
// SUMMARY : ...
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : Frederic Hecht
// E-MAIL  : frederic.hecht@sorbonne-universite.fr

#include "environment.hpp"
#include "iostream"
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <algorithm>

// set in getprog-unix.hpp in Graphic dir
const char *prognamearg = 0;
extern void (*initparallele)(int &, char **&); // to know if mpiversion ...

#ifdef PURE_WIN32
#include <windows.h>
#endif

#ifdef HAVE_GETENV
#include <cstdlib>
#endif
using namespace std;
bool load(string s);

const char SLACH = '/';
const char BACKSLACH = '\\';
#ifdef PURE_WIN32
const char dirsep = BACKSLACH, dirnsep = SLACH;
#else
const char dirnsep = BACKSLACH, dirsep = SLACH;
#endif

#include <sys/stat.h>

int dirExists (const string & path) {
  struct stat info;

  if (stat(path.c_str(), &info ) != 0)
    return 0;
  else if (info.st_mode & S_IFDIR)
    return 1;
  else
    return 0;
}

string DirName (const char *f) {
  const char *c = strrchr(f, dirsep);
  if (!c) return string("");
  else return string(f, strlen(f)-strlen(c));
}

string TransDir (string dir, string adddir="") {
  for (size_t i = 0; i < dir.size(); ++i)
    if (dir[i] == dirnsep) dir[i] = dirsep;
  if (dir.size() > 1 && dir[dir.size()-1] != dirsep)
    dir += dirsep;
  if (adddir.length() && dir[0] == '!')
    dir = adddir + dirsep + dir.substr(1);
  return dir;
}

template<typename T>
void show (const char *s, const T &l, const char *separateur="\n") {
  cout << s << * separateur;
  for (typename T::const_iterator i = l.begin(); i != l.end(); i++)
    cout << * i << * separateur;
}

bool EnvironmentFind (string key, string item) {
  EnvironmentData::iterator ekey = ffenvironment.find(key);
  if (ekey != ffenvironment.end()) {
    OneEnvironmentData *pl = &ekey->second;
    OneEnvironmentData::iterator i = find(pl->begin(), pl->end(), item);
    return i != pl->end();
  }

  return false;
}

bool EnvironmentClean (string key) {
  EnvironmentData::iterator ekey = ffenvironment.find(key);
  if (ekey != ffenvironment.end()) {
    OneEnvironmentData *pl = &ekey->second;
    pl->clear();
    return true;
  }

  return false;
}

bool EnvironmentInsert (string key, string item, string before) {
  bool ret = true;
  OneEnvironmentData &l = ffenvironment[key];
  char sufmpi[] = {'m', 'p', 'i', dirsep, '\0'};
  string suf = ((key == "loadpath") && initparallele ) ? sufmpi : "";
  if (verbosity > 1000) cout << " **EnvironmentInsert " << initparallele << " suf '"
    << suf << "' " << item << endl;
  if (!suf.empty() && dirExists(item+suf)) {
    if (verbosity >= 100) cout << " EnvironmentInsert: Add suf " << suf << " to " << item << " in GetEnvironment " << key << endl;
    item += suf;
  }

  OneEnvironmentData::iterator i = find(l.begin(), l.end(), item);

  if (i != l.end()) { ret = false; l.erase(i); } // if exist => remove
  i = find(l.begin(), l.end(), before);
  if (verbosity >= 100) cout << " insert " << key << " " << item << " " << before << endl;
  if (i == l.end() && before != "$")
    l.insert(l.begin(), item); // insert in front
  else
    l.insert(i, item); // insert before i

  return ret;
}

int GetEnvironment (const string &key, string items) {
  if (verbosity >= 100) cout << key << " -> " << items << endl;
  bool path = key.find("path") != string::npos;
  int d = 0, k = 0;
  if (path)
    items += ";;";
  for (size_t i = 0; i < items.size(); i++)
    if (items[i] == ';') {
      string item = items.substr(d, i-d);
      if (path) item = TransDir(item);
      if (verbosity >= 100) cout << " + " << item << endl;
      if (!EnvironmentFind(key, item)) {
        EnvironmentInsert(key, item, "$");
        k++;
      }
      d = i + 1;
    }

  return k;
}

int readinitfile (const string &file) {
	string line="";
	string key;
	string value;
	bool add=false;
	ifstream f(file.c_str());
	if( ! f) {
		if(verbosity>=20) cout << "error opening init  file: " << file << endl;
		return 0;}
        string dirfile=DirName(file.c_str());
	char c,c1,bv=0;
	int linenumber = 1;
	bool inkey=false,invalue=false,cmm=false;
	while (f)
	{
		c= f.get();
		c1=f.peek();
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
                {
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
					EnvironmentInsert(key,TransDir(value,dirfile),"$");
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
	}
	if( inkey || invalue || bv )
	  {
	    cout << " error read init file : " << file << " " <<  linenumber << endl;
	    cout << " line : " << line << endl;
	    return -11;
	  }
	return 1;
}


void GetEnvironment () {
 char  * ff_verbosity=0,* ff_loadpath=0,* ff_incpath=0,* home=0;

 // FFCS: we must make sure that FFCS does not reuse the same freefem++.pref as FF because some shared libraries must be
 // recompiled.
#ifdef PURE_WIN32
string     ffprefsuffix ="pref";
#else
string     ffprefsuffix ="pref";
#endif
#ifdef ENABLE_FFCS
 string  ffpref="freefem++-cs.pref";
#else
  string  ffpref="freefem++."+ffprefsuffix;
#endif

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
  char execpath[MAX_PATH+1];

  if (GetEnvironmentVariable("FF_VERBOSITY", envv, LEN) > 0)
   ff_verbosity=envv;

  if (GetEnvironmentVariable("FF_LOADPATH", envl, LEN) > 0)
   ff_loadpath=envl;

  if (GetEnvironmentVariable("FF_INCLUDEPATH", envi, LEN) > 0)
   ff_incpath=envi;

  if (GetEnvironmentVariable("HOMEPATH", envh, LEN) > 0)
    home=envh;
#endif
  if ( ff_verbosity ) {

    verbosity = atoi(ff_verbosity);
      if(verbosity > 4)
    cout << " -- GetEnvironmentVariable: verbosity= " <<verbosity   << endl;

  }

#ifdef PURE_WIN32
     int bytes = GetModuleFileName(NULL, execpath, MAX_PATH);
     execpath[bytes]='\0';
     if(bytes)
     {
         string execdir=DirName(execpath);
         if(execdir.length())
         {
             EnvironmentInsert("init-files",execdir+"\\"+ffpref,"$");
         }
     }
#else
  EnvironmentInsert("init-files","/etc/"+ffpref,"$");
#ifdef FF_PREFIX_DIR_APPLE
  EnvironmentInsert("init-files",string(FF_PREFIX_DIR_APPLE) + "/etc/" + ffpref ,"$");
#endif
#ifdef FF_PREFIX_DIR
  EnvironmentInsert("init-files",string(FF_PREFIX_DIR) + "/etc/" + ffpref  ,"$");
#endif

  if(prognamearg)
    {
    if( strchr(prognamearg,dirsep) )
      {
	EnvironmentInsert("init-files",TransDir(DirName(prognamearg))+"/../etc/"+ffpref,"$");
      }
    }
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

  if(ff_loadpath)
    GetEnvironment("loadpath",ff_loadpath);
  if(ff_incpath)
    GetEnvironment("includepath",ff_incpath);

   EnvironmentInsert("includepath","","");//    always add "" the first include path
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
    if(verbosity>10) cout << " --  GetEnvironment: verbosity is set to " << verbosity  << endl;

 }
const char *check_plugin=0;
void EnvironmentLoad() {
    if(check_plugin) {
        bool ok=load(check_plugin);
        if(ok) exit(0);
        else exit(1); }
    EnvironmentData::iterator toload=ffenvironment.find("load");
    if(  toload != ffenvironment.end())

	for (OneEnvironmentData::iterator i=toload->second.begin(); i != toload->second.end(); ++i)
	{
	    if(verbosity) cout << "PreEnv load :"<< *i << endl;
	    load(*i);
	}

}

// from ffapi to env. F. Hecht ..
#include "ffapi.hpp"
#include <dirent.h>
#include <strings.h>
#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif

namespace ffapi
{
// to change to tmp dir for exec ...
long chtmpdir()
{
    char tmp[256];
#ifdef _WIN32
    strcpy(tmp,"c:\\Temp");
    if (GetEnvironmentVariable("TEMP", tmp, 256) > 0);
#else
    strcpy(tmp,"/tmp/");
#endif
    if(verbosity>2)
        std::cout << " Change to " << endl;
    return chdir(tmp);

}
bool ff_justcompile=false;
bool ff_ch2edpdtmpir=0;
void ifchtmpdir()
{
    if(ff_ch2edpdtmpir) {
    }
}
}

#ifdef TESTMAIN
long verbosity = 50;
EnvironmentData environment;
int main () {
  GetEnvironment();
  return 0;
}
#endif
