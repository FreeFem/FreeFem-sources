//#pragma dont_inline on
//#pragma inline_depth(1)

#include "config-wrapper.h"

#include <complex>
#include "AFunction.hpp"
#include <cstdarg>
#include <cstring>
#include "error.hpp"
#include "lex.hpp"

#include "RNM.hpp"

#include "Operator.hpp"
// for exec routine 
#include "rgraph.hpp"
#include "InitFunct.hpp"

queue<pair<const E_Routine*,int> > debugstack;


double  VersionNumber(); 

OneOperator::pair_find OneOperator::Find(const ArrayOfaType & at)const
 { 
      const OneOperator *w=0,*oo;
      int nn=0,p=0;
 /*     for (oo=this;oo;oo=oo->next)
        if (oo->pref>=p && oo->WithOutCast(at)) 
          {
           if(p<oo->pref) {nn=0;p=oo->pref;}
           nn++;
           w=oo;}
      if (nn) return make_pair(w,nn);*/
      for (int ncast=0;ncast<=n;ncast++) // loop on the number of cast 
       {
         p=0;
         for (oo=this;oo;oo=oo->next)
          if (oo->pref>=p && oo->WithCast(at,ncast)) 
          { 
           if(p<oo->pref) {nn=0;p=oo->pref;}
            nn++;
            w=oo;}
         if (nn) return make_pair(w,nn);
       }
      for (oo=this;oo;oo=oo->next)
        if (oo->WithCast(at)) 
          {nn++;
           w=oo;}
       return make_pair(w,nn);       
}

OneOperator::pair_find OneOperator::FindWithOutCast(const ArrayOfaType & at)const
 { 
      const OneOperator *w=0,*oo;
      int n=0;
      for (oo=this;oo;oo=oo->next)
        if (oo->WithOutCast(at)) 
          {n++;
           w=oo;}
      return make_pair(w,n);
}

OneOperator* OneOperator::FindSameR(const ArrayOfaType & at)
 { 
     if (!this) return 0;
      OneOperator *w=0,*oo,*r;
      int n=0;
      for (oo=this;oo;oo=oo->next)
        //if (oo->WithOutCast(at)) 
        if  (at==*oo)  n++,r=oo;        
        else if (oo->WithOutCast(at)) n++,r=oo;
     // if (n>1) cout << "FindSameR " << n << endl;
      return n==1 ? r : 0;
}
      
void OneOperator::Show(ostream &f) const
{         
   const OneOperator *w=0,*oo;
   for (oo=this;oo;oo=oo->next)
     f << "\t (" <<  *oo << ")\n";
 }   

void OneOperator::Show(const ArrayOfaType & at,ostream &f) const
{         
         const OneOperator *w=0,*oo;
         int n=0,np=0;
         for (oo=this;oo;oo=oo->next)
           if (oo->WithOutCast(at)) {n++;f << "\t (" <<  *oo << ")\n";}
         if(n==0) 
          for (oo=this;oo;oo=oo->next)
           if (oo->WithCast(at)) {
              n++;
              if (oo->pref) np++;
              if (oo->pref) 
                f <<   " c(" << oo->pref << ") \t (" <<  *oo << ")\n" ;
                else f <<  " \t c(" <<  *oo << ")\n";
              }
         if (n==0) 
          {
           f << " List of choices "<< endl;
           Show(f);            
          }
         else if (np != 1) 
           f << " We have ambiguity " << n << endl; 
 }   
       
const  OneOperator * Polymorphic::Find(const char *op, const  ArrayOfaType &at) const
  {
    const_iterator i=m.find(op);
    if (i!=m.end())  
      { 
       OneOperator::pair_find r=i->second->Find(at);
       if (r.second==1) return r.first;
       }    
    return 0;
  }
const  OneOperator * Polymorphic::FindWithOutCast(const char *op, const  ArrayOfaType &at) const
  {
    const_iterator i=m.find(op);
    if (i!=m.end())  
      { 
       OneOperator::pair_find r=i->second->FindWithOutCast(at);
       if (r.second==1) return r.first;
       }    
    return 0;
  }
 
  
void Polymorphic::Show(const char *op,const ArrayOfaType & at,ostream &f)  const
    {
    const_iterator i=m.find(op);
    if (i==m.end()) f << " unknow " << op << " operator " << endl;
    else i->second->Show(at,f);
  }

C_F0::C_F0(const Polymorphic * poly,const char *op,const basicAC_F0 & p)
{
  ArrayOfaType at(p); 
  if (poly) { // a Polymorphic => polymorphisme
      const  OneOperator *  ff=poly->Find(op,at);
      if (ff) { 
        /* cout << endl;
         poly->Show(op,at,cout);
         cout << op << ": (in " << at << ") => " << " " << *ff<< "\n\n";*/
         f=ff->code(p);
         r=*ff;}
      else
        { 
           cerr << " error operator " << op << " " << at << endl;
           poly->Show(op,at,cerr);
           CompileError();
        }
   }
  else { 
      //  no polymorphisme
     cerr << " const Polymorphic * poly,const char *op,const basicAC_F0 & p)   " << endl;
     cerr  << op << " " << p << endl;
     CompileError();          
  }
}
       



//  operator without parameter
C_F0::C_F0(const Polymorphic * pop,const char *op) 
{
  basicAC_F0  p;
  p=0;
  *this= C_F0(pop,op,p);
}
//  operator unaire
C_F0::C_F0(const Polymorphic * pop,const char *op,const C_F0 & aa) 
{
  basicAC_F0  p;
  C_F0 a(aa);
  p=a;
  *this= C_F0(pop,op,p);
}
//  operator binaire
C_F0::C_F0(const Polymorphic * pop,const char *op,const  C_F0 & a,const  C_F0 & b) 
{
  C_F0 tab[2]={a,b};
  basicAC_F0  p;
  p=make_pair<int,C_F0*>(2,tab);
  *this= C_F0(pop,op,p);
}

//  operator trinaire
C_F0::C_F0(const Polymorphic * pop,const char *op,const  C_F0 & a,const  C_F0 & b,const  C_F0 & c) 
{
  C_F0 tab[3]={a,b,c};
  basicAC_F0  p;
  p=make_pair<int,C_F0*>(3,tab);
  *this= C_F0(pop,op,p);
}


 OneOperator::~OneOperator(){ 
       OneOperator * d=next;
       next=0; 
       if(! CodeAlloc::cleanning) // hash FH (pour les fuite de mémoire)
         while(d) 
        { 
         OneOperator * dd=d->next;
         d->next=0;
         delete d;
         d=dd;
        }
  }

void Polymorphic::Addp(const char * op,Value pp,...) const 
{
  pair<iterator,bool>  p=m.insert(make_pair<const Key,Value>(op,pp));
  Value f= p.first->second;	
  if (!p.second)  // not insert => old 
    *f += *pp;
  va_list ap;
  va_start(ap,pp);
  for(pp=va_arg(ap,OneOperator * );pp;pp=va_arg(ap,OneOperator * ))
    *f += *pp;
/*  if ( ! strlen(op) )  
   { // no polymorphisme
     if(m.size() !=1 ||  !f->Simple()) {
       cerr << " no polymorphisme and polymorphisme are mixed " << endl;
    //   for_each(m.begin,m.end(),ShowOn_cerr);
       CompileError();       
     }   
   } */
}

void Polymorphic::Add(const char * op,Value *pp) const
{
  if (*pp)
   {
    pair<iterator,bool>  p=m.insert(make_pair<const Key,Value>(op,*pp));
    Value f= p.first->second;	
    if (!p.second)  // not insert => old 
      *f += **pp;
    pp++;
    for(;*pp;pp++)
     *f += **pp;
   /*if ( ! strlen(op) )  
     { // no polymorphisme
      if(m.size() !=1 ||  !f->Simple()) {
       cerr << " no polymorphisme and polymorphisme are mixed " << endl;
      //   for_each(m.begin,m.end(),ShowOn_cerr);
       CompileError();       
     }  } */ }
     
}


 int  FindType(const char * name)  
   {
   C_F0 r;

     ListOfTOfId::const_iterator i=tables_of_identifier.begin();
      for(;i!=tables_of_identifier.end();++i)
      { 
      TableOfIdentifier * ti=*i;     
      r = ti->Find(name);
      if (r.NotNull()) return r.TYPEOFID();
    }
     return 0;
   } 

  C_F0 Find(const char * name)   
{
   C_F0 r;
   ListOfTOfId::const_iterator i=tables_of_identifier.begin();
   for(;i!=tables_of_identifier.end();++i)
    { 
      TableOfIdentifier * ti=*i;     
      r = ti->Find(name);
      if (r.NotNull()) return r;
    }
    cerr << " The Identifier " << name << " does not exist " << endl;
    CompileError();
    return r;
}

C_F0 TableOfIdentifier::destroy()
{
 int k=0;
// cout << "\n\t List of destroy variables " << m.size() << " : " ; 
 for (pKV * i=listofvar;i;i=i->second.next)
   {
     if  (i->second.del && i->second.first->ExistDestroy() )
    // cout  << i->first << ", " ;
     assert(i->second.first);
     if (i->second.del && i->second.first->ExistDestroy() ) k++;
   }
// cout << endl;
 ListOfInst *l=new ListOfInst(k);
 for (pKV * i=listofvar;i;i=i->second.next)
     if (i->second.del && i->second.first->ExistDestroy()) 
       l->Add(i->second.first->Destroy(i->second) );
   
 return C_F0(l);     
}


Expression basicForEachType::Destroy(const C_F0 & e) const 
{
    return destroy ? NewExpression(destroy,e) : (Expression)  e;
} 

basicForEachType::~basicForEachType()
  {
   if(casting) delete casting;
   ti.clear();
  }

 void basicForEachType::SetArgs(const ListOfId *lid) const
{ SHOWVERB(cout << "SetArgs::\n ") ;throwassert(lid==0 || lid->size()==0);}


const  Type_Expr &   TableOfIdentifier::New(Key k,const Type_Expr & v,bool del)
 {
    pair<iterator,bool>  p=m.insert(pKV(k,Value(v,listofvar,del)));
    listofvar = &*m.find(k);
    if (!p.second) 
     {
       cerr << " The identifier " << k << " exist \n";
       cerr << " \t  the existing type is " << *p.first->second.first << endl;
       cerr << " \t  the new  type is " << *v.first << endl;
       CompileError();
     }
    return v;
 }
 void  TableOfIdentifier::Add(Key k,Key op,OneOperator *p0,OneOperator *p1,
      OneOperator *p2,OneOperator *p3,OneOperator *p4,OneOperator *p5,OneOperator *p6)
{
  iterator i= m.find(k);
  if (i==m.end()) // new
    {
     Value poly0=Value(atype<Polymorphic*>(),new Polymorphic(),listofvar);     
    i=m.insert(make_pair<const Key,Value>(k,poly0)).first;
    listofvar= &*i;
    }
    const Polymorphic * p= dynamic_cast<const Polymorphic *>(i->second.second);
  if ( !p) { 
    cerr << k << " is not a Polymorphic id " << endl;
    CompileError();
  }
  p->Add(op,p0,p1,p2,p3,p4,p5,p6);
}

 ArrayOfaType::ArrayOfaType(const ListOfId * l)
  : n(l->size()),t(new aType[n]),ellipse(false)
 {
    for (int i=0;i<n;i++) 
      {
      t[i]=(*l)[i].r;  
       if ( ! t[i])
        {
           cerr << " Argument " << i << " '"<< (*l)[i].id << "' without type\n";
           CompileError("DCL routine: Argument without type ");
         }
      }
 }

bool ArrayOfaType::WithOutCast( const ArrayOfaType & a) const 
 {
   if ( ( !ellipse && (a.n != n))  || (ellipse && n > a.n) ) return false;
   for (int i=0;i<n;i++)
       if (! a.t[i]->SametypeRight(t[i])) 
        return false;       
 // cerr << " TRUE " << endl;    
   return true;
 }

  
bool ArrayOfaType::WithCast( const ArrayOfaType & a,int nbcast) const 
 {  
   if (  ( !ellipse && (a.n != n))  || (ellipse && n > a.n) ) return false;
   for (int i=0;i<n;i++)
     if ( a.t[i]->SametypeRight(t[i])) ;
     else if (! t[i]->CastingFrom(a.t[i])) return false; 
     else if ( --nbcast <0) return false;
   return true;
 }    
 
void basicForEachType::AddCast(CastFunc f1,CastFunc f2,CastFunc f3,CastFunc f4,
  CastFunc f5,CastFunc f6,CastFunc f7,CastFunc f8)
 {
   CastFunc ff[]={f1,f2,f3,f4,f5,f6,f7,f8,0};
   for (int i=0;ff[i];i++)
    {
      throwassert(this == *ff[i] );
      if (casting->FindSameR(*ff[i]))
       {
         cerr << " The casting to " << *ff[i] << " exist " << endl;
         cerr << " List of cast " << endl;
         casting->Show(cerr);
         CompileError();
       }
      if (casting)  *casting += *ff[i];
      else casting = ff[i];
/*      
   if( ! mapofcast.insert(make_pair<const aType,CastFunc>(ff[i]->a,ff[i])).second) 
     {
        cerr << " The casting to "<< *this << " from " << ff[i]->a << " exist " << endl;
        cerr << " List of cast " << endl;
        for_each(mapofcast.begin(),mapofcast.end(),CerrCast);
        CompileError();
     } */
  }
}

 ostream & operator<<(ostream & f,const OneOperator & a)     
{
//   for(const OneOperator * tt=&a;tt;tt=tt->next)
     f << "\t  " << * (a.r) << " :  "  <<(const ArrayOfaType &) a;
   return f;   
}

 ostream & operator<<(ostream & f,const Polymorphic & a)
{
  Polymorphic::const_iterator i;
  if(!&a) return f << "Null " << endl;
  for (i=a.m.begin();i!=a.m.end();i++)
   {
    f << "   operator" << i->first << " : " << endl;
    i->second->Show(f);
   }
  return f;
}  
 ostream & operator<<(ostream & f,const ArrayOfaType & a)
   { 
     for (int i=0;i<a.n;i++) 
       f <<  (i ? ", " : " ") << *a.t[i];
       if (a.ellipse ) f << "... ";
       else            f << " "; 
      return f;}
    ostream & operator<<(ostream & f,const TableOfIdentifier & t )
 {
   TableOfIdentifier::const_iterator i;
   for(i=t.m.begin();i!=t.m.end();i++)
    {
      TableOfIdentifier::Value v=i->second;
      f << i->first << ":  " << *v.first << " <- " ; 
      const Polymorphic * p=dynamic_cast<const Polymorphic *>(v.second);
      if(p) f << "Polymorphic " << *p << endl;
      else  f << " Simple @" <<  v.second << endl;
    }
    return f;
 }
 
Expression NewExpression(Function1 f,Expression a)
{
  throwassert(f);
  return new E_F0_Func1(f,a);
}
Expression NewExpression(Function2 f,Expression a,Expression b)
{
  throwassert(f);
  return new E_F0_Func2(f,a,b);
 
}
 
 void ShowType(ostream & f)
 {
 
   map<const string,basicForEachType *>::const_iterator i;
   for(i=map_type.begin();i!=map_type.end();i++)
     {
       f << " --"<< i->first <<" = " ; 
       i->second->Show(f) ;
       f << endl;
     }

 }

 void basicForEachType::Show(ostream & f) const {
       f << " " <<* this << endl;
       if (casting) casting->Show(f) ;
       if (ti.m.size())
        { 
          TableOfIdentifier::const_iterator mc=ti.m.begin();
          TableOfIdentifier::const_iterator end=ti.m.end();
          for (;mc != end;mc++)
          {
            f  << "    " << mc->first << ",  type :" <<  *mc->second.first << endl;
            const Polymorphic * op =dynamic_cast<const Polymorphic *>(mc->second.second) ;
            if ( op )  f << *op << endl;
          }
        }
   }
 
//size_t CodeAlloc::nb=0;
//size_t CodeAlloc::lg=0;
/*
FuncForEachType::FuncForEachType(const basicForEachType * t)
: rtype(t),basicForEachType(t->,sizeof((void*)()),
         0,0,0,0){ throwassert(t);}
*/

E_Routine::E_Routine(const Routine * routine,const basicAC_F0 & args)
 : rt(routine->tret),code(routine->ins),clean(routine->clean),
   nbparam(args.size()),param(new Expression[nbparam]),name(routine->name)
{    
   assert(routine->ins); 
   for (int i=0;i<args.size();i++)
        param[i]=routine->param[i].r->CastTo(args[i]);
};

AnyType E_Routine::operator()(Stack s)  const  {
   debugstack.push(make_pair<const E_Routine*,int>(this,TheCurrentLine));
   const int lgsave=BeginOffset*sizeof(void*);
   char  save[lgsave];
   AnyType ret=Nothing;
   memcpy(save,s,lgsave); // save register 
   AnyType *listparam=new AnyType[nbparam];
   for (int i=0;i<nbparam;i++)
     listparam[i]= (*param[i])(s); // set of the parameter 
   Stack_Ptr<AnyType>(s,ParamPtrOffset) = listparam; 
   try {  (*code)(s);  }
   catch( E_exception & e) { 
           (*clean)(s); 
          // cout << " catch " << e.what() << " clean & throw " << endl;
           if (e.type() == E_exception::e_return)  
              ret = e.r;
           else 
              ErrorExec("E_exception: break or contine not in loop ",1);
        }
        
   delete [] listparam; 
    memcpy(s,save,lgsave);  // restore register
    TheCurrentLine=debugstack.front().second;
    debugstack.pop();
 
   return ret;
}

void ShowDebugStack()
 {
   if (mpisize)
   cerr << "  current line = " << TheCurrentLine  
        << " mpirank " << mpirank << " / " << mpisize <<endl; 
   else
   cerr << "  current line = " << TheCurrentLine  << endl; 
   while ( debugstack.size() )
     {
        
        cerr << " call " << debugstack.front().first->name<< "  at  line " 
             <<debugstack.front().second << endl; 
        debugstack.pop();     
     }
 }
 

  int  E_F0::Optimize(deque<pair<Expression,int> > &l,MapOfE_F0 & m, size_t & n)  
     {
      int rr = find(m);
      if (rr) return rr;
      if( (verbosity / 10)% 10 == 1) 
      	 cout << "\n new expression : " << n  << " mi=" << MeshIndependent()<< " " << typeid(*this).name()  
      	      << " :" << *this << endl;        
       return insert(this,l,m,n);
     }   
     

class E_F0para :public E_F0 { public:
  const int i;
  AnyType operator()(Stack s)  const  {
  //  return  (* Stack_Ptr<Expression>(s,ParamPtrOffset)[i])(s);
    return Stack_Ptr<AnyType>(s,ParamPtrOffset)[i];
  }
   E_F0para(int ii) : i(ii){}
};        
     
Routine::Routine(aType tf,aType tr,const char * iden,  ListOfId *l,Block * & cb)
    : OneOperator(tr,l),offset(cb->OffSet(sizeof(void*))),
     tfunc(tf),tret(tr),name(iden),param(*l),
      currentblock(new Block(cb)),ins(0),clean(0) 
     {
       delete l;  // add  FH 24032005 (trap ) 
       cb = currentblock; 
       for (int i=0;i<param.size();i++)
           currentblock->NewID(param[i].r,param[i].id,C_F0(new E_F0para(i),param[i].r),!param[i].ref);
     }
   Block * Routine::Set(C_F0 instrs) 
       { 
         ins=instrs;
         clean = (C_F0) currentblock->close(currentblock);
         return    currentblock;} 
         
 
E_F0 * Routine::code(const basicAC_F0 & args) const 
{
   
   return new E_Routine(this,args);
}

void basicAC_F0::SetNameParam(int n,name_and_type *l , Expression * e) const 
{
 int k=0;
 if ( !n && !named_parameter)  return;
 
  for (int i=0;i<n;i++)
  {
     C_F0  ce=find(l[i].name) ;
     if (ce.LeftValue()==0)
       e[i]=0;
     else  {
       e[i]= map_type[l[i].type->name()]->CastTo(ce);
       k++;
       }
  }
  
 if (!named_parameter) return;
  
  if (k!=  named_parameter->size()) 
   {
      cout << " Sorry some name parameter are not used!  found" <<  k << " == " << named_parameter->size() <<endl;
      for(const_iterator ii=named_parameter->begin(); ii != named_parameter->end();ii++)
       {
        for (int i=0;i<n;i++)
          if (!strcmp(l[i].name,ii->first)) 
            goto L1;
         cout << "\t the parameter is '" << ii->first << "' is unused " << endl;
        L1:;
       }
    if ( n) {
    cerr << " The named parameter can be " << endl;
    for (int i=0;i<n;i++)
       cerr << "\t" << l[i].name << " =  <" << l[i].type->name() << ">\n";
    }
    CompileError("Unused named parameter");  
   }
}
