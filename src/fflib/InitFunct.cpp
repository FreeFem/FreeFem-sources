#include "InitFunct.hpp"
#include <algorithm>
#include <deque>
#include <iostream>

using namespace std;
typedef void  (* afunc)(); 
typedef pair<int,afunc> InitFunct;

deque<InitFunct> * getInitFunctlist()
{
  static deque<InitFunct> * data = new deque<InitFunct>();
  return data;
}


void call(const InitFunct & a) { (*a.second)();}
void show(const InitFunct & a) { cout << a.first << " " <<a.second << endl;}

bool comp(const InitFunct a,const InitFunct b)
 { 
   return a.first < b.first;
 }
 
 void  callInitsFunct() 
 {
   deque<InitFunct> *  l(getInitFunctlist()); 
   sort(l->begin(),l->end(),comp);
 //  cout << " callInitsFunct : " << l->size() << endl;
   for_each(l->begin(),l->end(),show);   
   for_each(l->begin(),l->end(),call);
 }
 
void  addInitFunct(int i,void  (* f)()) 
{ 
  getInitFunctlist()->push_back(make_pair(i,f));
//  cout << " addInitFunct: " << i << " " << f << endl; 
}
