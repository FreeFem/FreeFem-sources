#ifndef INITSFUNCT
#define INITSFUNCT_HPP_

void  addInitFunct(int i,void  (* f)()) ;
void  callInitsFunct() ;
struct  addingInitFunct { 
  addingInitFunct(int i,void  (* f)()) { addInitFunct(i,f);}
} ;
#endif
