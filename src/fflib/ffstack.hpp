//-----------------------------------  
//  to manage the freefem stack 

//   in the stack we save all the variable 
//   a adresse 0 we have the MeshPointStack to defineP,N, ....
//   a adresse sizeof(void *) 
// 
//
//  Offset in (void *)
const int MeshPointStackOffset =0;
const int ParamPtrOffset = 1;
const int ElemMatPtrOffset = 2;
const int ExprPtrs = 4;
const int NbPtrs = 5;



const int BeginOffset = 6;
//  0 : MeshPoint pointeur 
//  1 : ParamPtrOffset
//  2 : Truc les matrice elementaire

#define NEWFFSTACKxxx

#ifndef NEWFFSTACK

typedef void *Stack;


static Stack  NullStack=0;
//typedef StackType& Stack;


template<class T>
T * Stack_offset (Stack stack,size_t offset)  
  {return   (T *) (void *) (((char *) stack)+offset);}

template<class T>
T * & Stack_Ptr (Stack stack,size_t offset)  
  {return   (T * &)  (((void **) stack)[offset]);}

void ShowType(ostream & f);


struct VOIDPtrType {
  virtual ~VOIDPtrType();
};

template<class T> 
struct PtrType: public VOIDPtrType {
 T * p;
 PtrType(T* pp):p(pp) {}
 ~PtrType() {delete p;}
};

template<class T> 
struct PtrArrayType: public VOIDPtrType {
 T * p;
 PtrArrayType(T* pp):p(pp) {}
 PtrArrayType() {delete [] p;}
};


inline Stack newStack(size_t l)
 {
  Stack thestack = new char[l];
  for (int i = 0;i< l/sizeof(long);i++) ((long*) thestack)[i]=0;
  ((char **) thestack)[MeshPointStackOffset] = new char [1000]; 
  
  return thestack;
 // return *new StackType(l);}
}

inline void deleteStack(Stack s) 
 {
    delete [] (((char **)  s)[MeshPointStackOffset]);
    delete [] (char *) s;
 // s.clean();
 }
#else

struct StackType;

//typedef void *Stack;

typedef StackType & Stack;

struct StackType {
 size_t lg;
 char * stack; 
 char * MeshPointStack;
 operator char *() { return stack;}
 operator void *() { return stack;}
 operator long *() { return (long *)(void *)stack;}
 operator void **() {return (void **) (void *) stack;}
 template<class T> 
 T * Offset(size_t offset){ return (T*) (void*) (stack+offset);}
 template<class T> 
 T *& ptr(size_t offset){ return (T* &) ((void**) (void *) stack)[offset];}
 StackType(size_t ll) :lg(ll),stack(new char[ll]),MeshPointStack(new char[1000]) 
  {
  long * p= ptr<long>(0);
  long l4=lg/sizeof(long);
  for (int i = 0;i< l4;i++) p[i]=0;

  ptr<char>(MeshPointStackOffset)=MeshPointStack;
  }
 void clean() { delete []stack; delete [] MeshPointStack; }
};

StackType * NullStackPtr= 0;
StackType & NullStack(*NullStackPtr);
//typedef StackType& Stack;


template<class T>
T * Stack_offset (Stack stack,size_t offset)  
  {return   stack.Offset<T>(offset);}

template<class T>
T * & Stack_Ptr (Stack stack,size_t offset)  
  {return   (T * &)  (((void **) (void *) stack.stack)[offset]);}

void ShowType(ostream & f);


struct VOIDPtrType {
  virtual ~VOIDPtrType();
};

template<class T> 
struct PtrType: public VOIDPtrType {
 T * p;
 PtrType(T* pp):p(pp) {}
 ~PtrType() {delete p;}
};

template<class T> 
struct PtrArrayType: public VOIDPtrType {
 T * p;
 PtrArrayType(T* pp):p(pp) {}
 PtrArrayType() {delete [] p;}
};


inline Stack newStack(size_t l)
 {
/*  Stack thestack = new char[l];
  for (int i = 0;i< l/sizeof(long);i++) ((long*) thestack)[i]=0;
  ((char **) thestack)[MeshPointStackOffset] = new char [1000]; 
  
  return thestack;*/
  return *new StackType(l);
  }


inline void deleteStack(Stack s) 
 {
   // delete [] (((char **)  s)[MeshPointStackOffset]);
  //  delete [] (char *) s;
 s.clean();
 }
#endif
//------------------------------------
