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


 const Stack  NullStack=0;
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



#endif
//------------------------------------
 
 // Add FH mars 2006 
 // clean pointeur ceated by the language 
 // ----
 struct BaseNewInStack {
  virtual ~BaseNewInStack() {};
};

 struct BaseNewInStack ;
 struct StackOfPtr2Free;

inline StackOfPtr2Free  * & WhereStackOfPtr2Free(Stack s) { return  Stack_Ptr<StackOfPtr2Free>(s,ExprPtrs) ;} // fait  


struct StackOfPtr2Free {
	typedef vector<BaseNewInStack *>::iterator iterator;
	StackOfPtr2Free  ** where; // where is store the ptr to the stack 
	StackOfPtr2Free *prev; // previous stack 

	vector<BaseNewInStack *> stackptr;
	void add(BaseNewInStack *p) { 
	   // cout << "\n\t\t ### ptr/lg add  " << p << "  at " << stackptr.size() << "  \n";	   
	   stackptr.push_back(p);}
	
public: 
	StackOfPtr2Free(Stack s):
		where(&WhereStackOfPtr2Free(s)),
		prev(*where)
	       { 
			stackptr.reserve(20); 	
			if(prev) Add2StackOfPtr2Free(s,this);
	       }
	
	bool clean() 
	 { bool ret= !stackptr.empty();
	    if(ret)
	      { 
	        if(stackptr.size()>=20 && verbosity>2) 
	           cout << "\n\t\t ### big?? ptr/lg clean " << stackptr.size() << " ptr's\n ";
		for (iterator i=stackptr.end(); i != stackptr.begin();)
		{
		   
			delete  (* (--i) ); 
		     //cout << "StackOfPtr2Free: clean " << (* (i) ) << endl;
		}
		stackptr.resize(0);// clean the
		
	     }
	   return ret;
	}
	
	~StackOfPtr2Free() {clean(); *where=prev;} // restore the previous stack
private:// no copy ....
 	StackOfPtr2Free(const StackOfPtr2Free&);
	void operator =(const StackOfPtr2Free&);

 template<class T>
  friend  T * Add2StackOfPtr2Free(Stack s,T * p);
	
};	

 


template<class T>
struct NewInStack: public BaseNewInStack   {	
   T * p;
   bool array;
  ~NewInStack() { 
     // cout << "~NewInStack " << typeid(T).name() << " " << p << "  array " << array << " " << this << "\n";
      if(p) 
        delete p;}  
private: 
   NewInStack(T * pp,bool aa=false) : p(pp),array(aa) {
      //  cout << "NewInStack " << typeid(T).name() << " "  << p << "  array " << aa << " " << this << "\n";
   } 
   
   
 template<class TT> 
 friend  TT * Add2StackOfPtr2FreeA(Stack s,TT * p);
 template<class TT> 
 friend  TT * Add2StackOfPtr2Free(Stack s,TT * p);
   
};
// ajout of 2 class NewRefCountInStack and NewArrayInStack
//  for clean of meshes 

template<class T>
struct NewRefCountInStack: public BaseNewInStack   {	
    T * p;
    bool array;
    ~NewRefCountInStack() { 
	// cout << "~NewInStack " << typeid(T).name() << " " << p << "  array " << array << " " << this << "\n";
	if(p) p->destroy();}  
private: 
    NewRefCountInStack(T * pp,bool aa=false) : p(pp),array(aa) {
	//  cout << "NewInStack " << typeid(T).name() << " "  << p << "  array " << aa << " " << this << "\n";
    } 
    
    template<class TT> 
    friend  TT * Add2StackOfPtr2FreeRC(Stack s,TT * p);
    
};

template<class T>
struct NewArrayInStack: public BaseNewInStack   {	
    T * p;
    bool array;
    ~NewArrayInStack() { 
	// cout << "~NewInStack " << typeid(T).name() << " " << p << "  array " << array << " " << this << "\n";
	if(p)  delete [] p;
   }  
private: 
    NewArrayInStack(T * pp,bool aa=false) : p(pp),array(aa) {
	//  cout << "NewInStack " << typeid(T).name() << " "  << p << "  array " << aa << " " << this << "\n";
    } 
    
    
    template<class TT> 
    friend  TT * Add2StackOfPtr2FreeA(Stack s,TT * p);
    template<class TT> 
    friend  TT * Add2StackOfPtr2Free(Stack s,TT * p);
    
};

template<class T>
T * Add2StackOfPtr2FreeRC(Stack s,T * p)
{
    if(p)	
	WhereStackOfPtr2Free(s)->add(new NewRefCountInStack<T>(p));
    return p;
}	

template<class T>
T * Add2StackOfPtr2Free(Stack s,T * p)
{
   if(p)	
     WhereStackOfPtr2Free(s)->add(new NewInStack<T>(p));
   return p;
}	
template<class T>
T * Add2StackOfPtr2FreeA(Stack s,T * p)
{
   if(p)	
     WhereStackOfPtr2Free(s)->add(new NewArrayInStack<T>(p));
   return p;
}	
//  fin modif gestion of allocation of Ptr in Language 
//  ---------------------------------------------------
#ifndef NEWFFSTACK

inline Stack newStack(size_t l)
 {
  Stack thestack = new char[l];
  for (size_t i = 0;i< l/sizeof(long);i++) ((long*) thestack)[i]=0;
  ((char **) thestack)[MeshPointStackOffset] = new char [1000]; 
  WhereStackOfPtr2Free(thestack)=new StackOfPtr2Free(thestack); 
  
  return thestack;
 // return *new StackType(l);}
}

inline void deleteStack(Stack s) 
 {
    delete  WhereStackOfPtr2Free(s); // add gestion of the Ptr
    delete [] (((char **)  s)[MeshPointStackOffset]);
    delete [] (char *) s;
 // s.clean();
 }
#else
  a faire ....
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
