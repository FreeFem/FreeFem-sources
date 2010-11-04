template<class T> struct  Resize{ T *v;
  Resize( T * vv) : v(vv) {}
 }; 

template<class T> T *resize1(const Resize<T> & t,const long &n)
 {  
  t.v->resize(n);
  return  t.v;
 }

template<class T> T *resizeandclean1(const Resize<T> & t,const long &n)
 {  
  int nn= t.v->N(); // old size 
  
  for (int i=n;i<nn;i++)  {delete (*t.v)[i];} // clean
  t.v->resize(n);
  for (int i=nn;i<n;i++)  {(*t.v)[i]=0;}  
  return  t.v;
 }
 
template<class T> T *resize2(const Resize<T> & t,const long &n, const long & m)
 {  
  t.v->resize(n,m);
  return  t.v;
 }

template<class T> Resize<T> to_Resize( T *v){ return Resize<T>(v);}
