#ifndef LABEL_HPP
#define LABEL_HPP
 
class Label {  // reference number for the physics
  friend inline ostream& operator <<(ostream& f,const Label & r  )
    { f <<  r.lab ; return f; }
  friend inline istream& operator >>(istream& f, Label & r  )
    { f >>  r.lab ; return f; }
  public: 
  int lab;
  Label(int r=0):lab(r){}
  bool onGamma() const { return lab;} 
  int operator!() const{return !lab;} 
  bool operator<(const Label & r) const {return lab < r.lab;} 
  bool operator==(const Label & r) const {return lab == r.lab;} 
  bool operator>(const Label & r) const { return lab > r.lab;} 

  };
#endif
