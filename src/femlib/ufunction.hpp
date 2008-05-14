#ifndef UFUNCTION_HPP_
#define UFUNCTION_HPP_
// some usefull function
template<class T> inline T Min (const T &a,const T &b){return a < b ? a : b;}
template<class T> inline T Max (const T &a,const T & b){return a > b ? a : b;}
template<class T> inline T Abs (const T &a){return a <0 ? -a : a;}
template<class T> inline double Norme (const T &a){return sqrt(a*a);}
template<class T> inline void Exchange (T& a,T& b) {T c=a;a=b;b=c;}
template<class T> inline T Max (const T &a,const T & b,const T & c){return Max(Max(a,b),c);}
template<class T> inline T Min (const T &a,const T & b,const T & c){return Min(Min(a,b),c);}


template<int d,typename A,typename B>
struct Copy {
static A * copy(A a[d],B * const b)
{
  for (int i=0;i<d;++i)
    a[i]=b[i];
  return a;
}};

template<typename A,typename B>
struct Copy<1,A,B> {
static A * copy(A a[1],B * const b  )
{
    a[0]=b[0];
  return a;
}};

template<typename A,typename B>
struct Copy<2,A,B> {
static A * copy(A a[2],B * const b)
{
    a[0]=b[0];
    a[1]=b[1];

  return a;
}};

template<typename A,typename B>
struct Copy<3,A,B> {
static A * copy(A a[3],B * const b)
{
    a[0]=b[0];
    a[1]=b[1];
    a[2]=b[2];
  return a;
}};

template<typename A,typename B>
struct Copy<4,A,B> {
static A * copy(A a[4],B * const b)
{
    a[0]=b[0];
    a[1]=b[1];
    a[2]=b[2];
    a[3]=b[3];
  return a;
}};

struct UniqueId {
private:
  static int count ;
  int id;
public:
  UniqueId() : id(++count) {}
  bool operator==(UniqueId u) const {return id==u.id;}
  bool operator!=(UniqueId u) const {return id!=u.id;}
  void init(){id=++count;}
};

#endif
