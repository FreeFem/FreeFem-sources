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
#include "InitFunct.hpp"
#include <algorithm>
#include <deque>
#include <iostream>
#include <set> 
#include "ffapi.hpp"  
using namespace std;
typedef void  (* afunc)(); 
typedef pair<int,afunc> InitFunct;

deque<InitFunct> * getInitFunctlist()
{
  static deque<InitFunct> * data = new deque<InitFunct>();
  return data;
}

extern long verbosity;
set<string> & ff_SetofInitFunct() { static set<string> sset; return sset;}
void call(const InitFunct & a) { 
  if(verbosity>5) 
    cout << "\n addInitFunct : " << a.first << " call : " <<a.second  << " ( " ; 
  (*a.second)();  
  if(verbosity>5)
    cout <<  " ) " ;
}
bool comp(const InitFunct a,const InitFunct b)
 { 
   return a.first < b.first;
 }
 
 void  callInitsFunct() 
 {
   deque<InitFunct> *  l(getInitFunctlist()); 
   sort(l->begin(),l->end(),comp);
   if(verbosity>5) cout << " callInitsFunct : " << l->size() << endl;
   //   for_each(l->begin(),l->end(),show);   
   for_each(l->begin(),l->end(),call);
    l->clear();
 }
 
void  addInitFunct(int i,void  (* f)(),const char *name) 
{

  if(!name || (! *name ) ||  ff_SetofInitFunct().insert(name).second)
    { 
    getInitFunctlist()->push_back(make_pair(i,f));
      cout << " -- addInitFunct: " << i << " " << f 
			  << " " <<  (name ? name : " -- " ) <<endl; 
    }
  else 
    cout << " ********  addInitFunct "<< name << " is always load (skip) !" << endl; 

}



