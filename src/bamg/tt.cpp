#include <new.h>
#include <iostream.h> 
class A{ public:
  int NbRef;
  int allocated;
  static int IsCreateByNew;
  A() :NbRef(0),allocated(IsCreateByNew){IsCreateByNew=0;cout << " constructed  " << NbRef << " " << allocated << endl;}
 void  AddRef() {NbRef++;}
  ~A(){IsCreateByNew=allocated;}
  void * operator new(size_t l){IsCreateByNew=1;cout << "new " << l  << endl; return ::operator new(l);}
  void operator delete(void *p,size_t l){cout << "delete " << p << " " << l<< endl;
  if(IsCreateByNew)::operator delete (p);
  else cerr << " Try to remove no allocated pointeur " <<p <<  endl;IsCreateByNew=0;}
};
 int A::IsCreateByNew=0;
int main()
{
  A a;
  A *b( new A()); 
  delete b;
  delete &a;
  return 0;
}
