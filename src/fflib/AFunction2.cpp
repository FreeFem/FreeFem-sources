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
//#pragma dont_inline on
//#pragma inline_depth(1)

// TODO: remove this block as soon as autoconf is removed from FreeFem++
#ifndef CMAKE
#include <config.h>
#endif

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

vector<pair<const E_Routine*,int> > *debugstack=0;


class vectorOfInst : public  E_F0mps { public:
    int n;
    Expression * v;
    vectorOfInst(int k): n(k),v(new Expression[k]) {ffassert(v);
      for(int i=0;i<n;++i) v[i]=0; }
    ~vectorOfInst(){ delete [] v;}
    bool empty() const {return n;}

   AnyType operator()(Stack s)  const {
     for (int i=0;i<n;++i)
      {
       ffassert(v[i]);
       (*(v[i]))(s);
      }
      return Nothing;
   }
   void  eval(Stack s,int j=-1)  const {
       if(j>=0) j = n-j;
       else j =0;
       if(verbosity>999) cout << " eval vectorOfInst " << j << " " << n << endl;
        for (int i=j;i<n;++i)
        {
            ffassert(v[i]);
            (*(v[i]))(s);
        }
    }

  private:
  vectorOfInst(const vectorOfInst &);
  void operator=(const vectorOfInst &);
};

double  VersionNumber();

OneOperator::pair_find OneOperator::Find(const ArrayOfaType & at)const
 {
      const OneOperator *w=0,*oo;
      int nn=0,p=-10000;
      for (int ncast=0;ncast<=n;ncast++) // loop on the number of cast
       {
         p=-10000;
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

// <<FindSameR>>
OneOperator* OneOperator::FindSameR(const ArrayOfaType & at)
 {
     if (this==tnull) return 0;
      OneOperator *oo,*r;
      int n=0;
      for (oo=this;oo;oo=oo->next)
        {
        if  (at==*oo)  n++,r=oo;
        else if (oo->WithOutCast(at)) n++,r=oo;
        }
      return n==1 ? r : 0;
}

void OneOperator::Show(ostream &f) const
{
   const OneOperator *oo;
   for (oo=this;oo;oo=oo->next)
     f << "\t (" <<  *oo << ")\n";
 }

void OneOperator::Show(const ArrayOfaType & at,ostream &f) const
{
         const OneOperator *oo;
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
      int nf=0;
    if (i!=m.end())
      {
       OneOperator::pair_find r=i->second->Find(at);
       if (r.second==1) return r.first;
          nf=max(nf,r.second);
       }
      if(nf) { cerr << "\n Warning ambiguity Polymorphic Find "<<  nf << endl;
          Show(op,at,cerr); }
    return 0;
  }
const  OneOperator * Polymorphic::FindWithOutCast(const char *op, const  ArrayOfaType &at) const
  {
      int nf=0;
    const_iterator i=m.find(op);
    if (i!=m.end())
      {
       OneOperator::pair_find r=i->second->FindWithOutCast(at);
       if (r.second==1) return r.first;
           nf=max(nf,r.second);
       }
      if(nf) { cerr << "\n Warning ambiguity Polymorphic FindWithOutCast "<<op<< " "<<  nf << endl;  Show(op,at,cerr);}
    return 0;
  }


void Polymorphic::Show(const char *op,const ArrayOfaType & at,ostream &f)  const
    {
    const_iterator i=m.find(op);
    if (i==m.end()) f << " unknow " << op << " operator " << endl;
    else i->second->Show(at,f);
  }

// <<C_F0_constructor_pop_char_basicAC_F0_impl>> cf [[file:AFunction.hpp::C_F0_constructor_pop_char_basicAC_F0_decl]]
C_F0::C_F0(const Polymorphic * poly,const char *op,const basicAC_F0 & p)
{
    ArrayOfaType at(p);
    if (poly) { // a Polymorphic => polymorphisme
	const  OneOperator *  ff=poly->Find(op,at);
	if (ff) {
            if( verbosity>9999) {cout << endl;
	     poly->Show(op,at,cout);
                cout << op << ": (in " << at << ") => " << " " << *ff<< "\n\n";}

	  // [[file:AFunction.hpp::OneOperator_code2]]
	  *this= ff->code2(p);
	}
	else
	  { if(mpirank==0)
	    {
		cerr << " error operator " << op << " " << at << endl;
		poly->Show(op,at,cerr);
		poly->Find(op,at);
	    }
	      CompileError();
	  }
    }
    else {
	//  no polymorphisme
	if(mpirank==0){
	    cerr << " const Polymorphic * poly,const char *op,const basicAC_F0 & p)   " << endl;
	    cerr  << op << " " << at << endl;
            
	}
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

// <<C_F0_constructor_binary_operator>> operator binaire
C_F0::C_F0(const Polymorphic * pop,const char *op,const  C_F0 & a,const  C_F0 & b)
{
  C_F0 tab[2]={a,b};
  basicAC_F0 p;
  p=make_pair<int,C_F0*>(2,tab);

  // [[file:AFunction.hpp::C_F0_constructor_pop_char_basicAC_F0_decl]]
  *this=C_F0(pop,op,p);
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
       if(! CodeAlloc::cleanning) // hash FH (pour les fuite de m�moire)
         while(d)
        {
         OneOperator * dd=d->next;
         d->next=0;
         delete d;
         d=dd;
        }
  }

    OneOperator::OneOperator(aType rr)
      : ArrayOfaType(),r(rr),next(0),pref(0) {throwassert(r);}
    OneOperator::OneOperator(aType rr,aType  a)
      : ArrayOfaType(a,false),r(rr),next(0),pref(0) {throwassert(rr && a );}
    OneOperator::OneOperator(aType rr,aType  a,aType  b)
      : ArrayOfaType(a,b,false),r(rr),next(0),pref(0) {
     throwassert(rr && a && b);}
    OneOperator::OneOperator(aType rr,aType  a,aType  b,aType c)
      : ArrayOfaType(a,b,c,false),r(rr),next(0),pref(0)
        {throwassert(rr && a && b && c);}
    OneOperator::OneOperator(aType rr,aType  a,aType  b,aType c,aType d)
      : ArrayOfaType(a,b,c,d,false),r(rr),next(0),pref(0)
      {throwassert(rr && a && b && c);}

    OneOperator::OneOperator(aType rr,aType  a,aType  b,aType c,aType d,aType e)
      : ArrayOfaType(a,b,c,d,e,false),r(rr),next(0),pref(0)
       {throwassert(rr && a && b && c && d);} // Added by Fabian Dortu (5 parameters)
    OneOperator::OneOperator(aType rr,aType  a,aType  b,aType c,aType d,aType e,aType f)
      : ArrayOfaType(a,b,c,d,e,f,false),r(rr),next(0),pref(0)
      {throwassert(rr && a && b && c && d && e && f);} // Added by Fabian Dortu (6 parameters)
    OneOperator::OneOperator(aType rr,aType  a,aType  b,aType c,aType d,aType e,aType f, aType g)
      : ArrayOfaType(a,b,c,d,e,f,g,false),r(rr),next(0),pref(0)
       {throwassert(rr && a && b && c && d && e && f && g);} // Added by Fabian Dortu (7 parameters)
    OneOperator::OneOperator(aType rr,aType  a,aType  b,aType c,aType d,aType e,aType f, aType g, aType h)
     : ArrayOfaType(a,b,c,d,e,f,g,h,false),r(rr),next(0),pref(0)
       {throwassert(rr && a && b && c && d && e && f && g && h);} // Added by Fabian Dortu (8 parameters)
    OneOperator::OneOperator(aType rr,aType  a,aType  b,aType c,aType d,aType e,aType f, aType g, aType h, aType i)
      : ArrayOfaType(a,b,c,d,e,f,g,h,i,false),r(rr),next(0),pref(0)
      {throwassert(rr && a && b && c && d && e && f && g && h && i);} // Added by Fabian Dortu (9 parameters)
    OneOperator::OneOperator(aType rr,aType  a,aType  b,aType c,aType d,aType e,aType f, aType g, aType h, aType i, aType j)
      : ArrayOfaType(a,b,c,d,e,f,g,h,i,j,false),r(rr),next(0),pref(0)
     {throwassert(rr && a && b && c && d && e && f && g && h && i && j);} // Added by Fabian Dortu (10 parameters)



    OneOperator::OneOperator(aType rr,const ArrayOfaType &ta)
      : ArrayOfaType(ta),r(rr),next(0),pref(0)
       {throwassert(rr);}
    OneOperator::OneOperator(aType rr,bool ellipse)
      : ArrayOfaType(ellipse),r(rr),next(0),pref(0)
        {throwassert(rr );}
    OneOperator::OneOperator(aType rr,const ListOfId *l)
      : ArrayOfaType(l),r(rr),next(0),pref(0)
      {throwassert(rr );}

void Polymorphic::Addp (const char *op, Value pp, ...) const {
  pair<iterator,bool> p = m.insert(pair<const Key, Value>(op, pp));
  Value f = p.first->second;
  if (!p.second)  // not insert => old
    *f += *pp;
  va_list ap;
  va_start(ap, pp);
  for (pp = va_arg(ap, OneOperator*); pp; pp = va_arg(ap, OneOperator*))
    *f += *pp;
  va_end(ap);
}

void Polymorphic::Add(const char * op,Value *pp) const
{
  if (*pp)
   {
    pair<iterator,bool>  p=m.insert(pair<const Key,Value>(op,*pp));
    Value f= p.first->second;
    if (!p.second)  // not insert => old
      *f += **pp;
    pp++;
    for(;*pp;pp++)
     *f += **pp;
 }

}


// <<FindType>>
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

/// <<Find>> uses [[file:global.cpp::tables_of_identifier]]

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
    if(mpirank==0)
    cerr << " The Identifier " << name << " does not exist " << endl;
    CompileError();
    return r;
}

vectorOfInst* TableOfIdentifier::newdestroy()
{
 int k=0;
 for (pKV * i=listofvar;i;i=i->second.next)
   {
     if  (i->second.del && i->second.first->ExistDestroy() )
     assert(i->second.first);
     if (i->second.del && i->second.first->ExistDestroy() ) k++;
   }
  ffassert(nIdWithDelete==k);
// new code
  vectorOfInst * l= new vectorOfInst(k);
  int j=0;
 for (pKV * i=listofvar;i;i=i->second.next)
     if (i->second.del && i->second.first->ExistDestroy())
       l->v[j++]=i->second.first->Destroy(i->second) ;
  ffassert(j==k);
 return l;
}
C_F0  TableOfIdentifier::destroy() {return C_F0(newdestroy());}
   void TableOfIdentifier::clear()
   {
     for (iterator i=m.begin();i!=m.end();++i)
       {
   //     delete i->first;
        }
     m.clear();
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

basicForEachType::basicForEachType(const type_info  & k,
                                          const size_t s,
                                          const E_F1_funcT_Type * p,
                                          basicForEachType *rr,
                                          Function1 iv,Function1 id, Function1  dreturn)
      : ktype(&k),//ktypefunc(0),
        size(s),
        un_ptr_type(rr?rr:this),
        casting(0), // no casting to
        un_ptr(p),
        InitExp(iv),
        DoOnReturn(dreturn),
        //funct_type(0),
        destroy(id) {}
 void basicForEachType::SetArgs(const ListOfId *lid) const
{ SHOWVERB(cout << "SetArgs::\n ") ;ffassert(lid==0 || lid->size()==0);}



 TableOfIdentifier::TableOfIdentifier() : listofvar(0),nIdWithDelete(0) {}
 TableOfIdentifier:: ~TableOfIdentifier() {}


Block::Block(Block * f):fatherblock(f),top(f?f->top:BeginOffset*sizeof(void*)),topmax(top)
    {
      itabl=tables_of_identifier.insert(tables_of_identifier.begin(),&table);
    }
Block::~Block(){}

   vectorOfInst * Block::snewclose(Block *& c) {
    Block * a=c;
    tables_of_identifier.erase(a->itabl);
    c=a->fatherblock;
    if (a->fatherblock) {a->fatherblock->topmax=a->topmax;
        a->fatherblock->top=a->top;}

    vectorOfInst * r;
    r = a->table.newdestroy();
    delete a;
    return r;}

CC_F0  Block::close(Block *& c,C_F0  ins)
{
    CListOfInst inst;
    CC_F0 cins; cins=ins;
    inst = cins;
    inst.setclose(Block::snewclose(c));
    CC_F0 rr;
    rr=inst;
    return rr;
}

 CC_F0  Block::close(Block *& c,CListOfInst  inst) {

     inst.setclose(Block::snewclose(c));
     CC_F0 rr;
     rr=inst;
     return rr;
 }

   Block * Block::open(Block *& cb)
   {
    Block *  ncb = new Block(cb);
    if(verbosity>99) cout << " Block::open  " << ncb <<  " " << cb << endl;

    return cb = ncb;
   }


const  Type_Expr &   TableOfIdentifier::New(Key k,const Type_Expr & v,bool del)
  {
    if( this != &Global) {
	if ( Global.m.find(k) != Global.m.end() )
	  {
	    if(mpirank==0 && (verbosity>0))
	      cerr << "\n *** Warning  The identifier " << k << " hide a Global identifier  \n";

	  }
    }
      pair<iterator,bool>  p=m.insert(pKV(k,Value(v,listofvar,del)));
      listofvar = &*m.find(k);
      if (!p.second)
	{
	    if(mpirank==0) {
		cerr << " The identifier " << k << " exists \n";
		cerr << " \t  the existing type is " << *p.first->second.first << endl;
		cerr << " \t  the new  type is " << *v.first << endl;
	    }
	    CompileError();
	}
      if(del && v.first->ExistDestroy() )
      {
          nIdWithDelete++;
          if(verbosity>9999) cout << "\n \t add ExistDestroy" << endl;

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
	    i=m.insert(pair<const Key,Value>(k,poly0)).first;
	    listofvar= &*i;
	}
      const Polymorphic * p= dynamic_cast<const Polymorphic *>(i->second.second);
      if ( !p) {
	  if(mpirank==0)
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
	   if(mpirank==0)
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
	    ffassert(this == *ff[i] );
	    if (casting->FindSameR(*ff[i]))
	      {
		  if(mpirank==0)
		    {
			cerr << " The casting to " << *ff[i] << " exists " << endl;
			cerr << " List of cast " << endl;
			casting->Show(cerr);
		    }
		  CompileError();
	      }
	    if (casting)  *casting += *ff[i];
	    else casting = ff[i];
	}
  }

 ostream & operator<<(ostream & f,const OneOperator & a)
{
     f << "\t  " << * (a.r) << " :  "  <<(const ArrayOfaType &) a;
   return f;
}

 ostream & operator<<(ostream & f,const Polymorphic & a)
{
  Polymorphic::const_iterator i;
    if(&a==E_F0::tnull) return f << "Null " << endl;
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
  ffassert(f);
  return new E_F0_Func1(f,a);
}
Expression NewExpression(Function2 f,Expression a,Expression b)
{
  ffassert(f);
  return new E_F0_Func2(f,a,b);

}

// <<ShowType>>
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





E_Routine::E_Routine(const Routine * routine,const basicAC_F0 & args)
  :    code(routine->ins),
       rt(routine->tret),
       nbparam(args.size()),
       param(new Expression[nbparam]),
       name(routine->name)
{
   assert(routine->ins);
   for (int i=0;i<args.size();i++)  //  bug pb copie des string   dec 2007  FH  ???????????????
   {
      if(verbosity>10000) cout << "E_Routine " << *routine->param[i].r << " <- " << *args[i].left() << endl;
        param[i]=routine->param[i].r->CastTo(args[i]);
   }
};

E_Routine::~E_Routine() {
    if(verbosity>10000) cout << "~E_Routine()"<< endl;
    delete [] param;}
AnyType E_Routine::operator()(Stack s)  const  {
   debugstack->push_back(pair<const E_Routine*,int>(this,TheCurrentLine));
   const int lgsave=BeginOffset*sizeof(void*);
   char  save[lgsave];
   AnyType ret=Nothing;
   memcpy(save,s,lgsave); // save
    AnyType *listparam;
    Add2StackOfPtr2FreeA(s,listparam=new AnyType[nbparam]);

//  to day the memory gestion of the local variable are static,
   for (int i=0;i<nbparam;i++)
     listparam[i]= (*param[i])(s); // set of the parameter
   Stack_Ptr<AnyType>(s,ParamPtrOffset) = listparam;
   WhereStackOfPtr2Free(s)=new StackOfPtr2Free(s);// FH mars 2006

   try {
      ret=(*code)(s);
     }
   catch( E_exception & e) {
            if (e.type() == E_exception::e_return)
              ret = e.r;
           else
              ErrorExec("E_exception: break or contine not in loop ",1);
  }
  catch(...) { // clean and rethrow the exception
      WhereStackOfPtr2Free(s)->clean(); // FH mars 2005
      memcpy(s,save,lgsave);  // restore
      TheCurrentLine=debugstack->back().second;
      debugstack->pop_back();

      throw ;
     }

  //  (*clean)(s); //  the clean is done in CleanE_Routine delete .
   //  delete [] listparam; after return
    memcpy(s,save,lgsave);  // restore
    TheCurrentLine=debugstack->back().second;
    debugstack->pop_back();

   // il faudrait que les variable locale soit detruire apres le return
   // cf routine clean, pour le cas ou l'on retourne un tableau local.
   // plus safe ?????  FH.  (fait 2008)
   // mais pb si   a = f()+g()   OK les pointeurs des instruction sont detruit
    //  en fin d'instruction programme de l'appelant  FH 2007
   // ... ou alors changer le return ???? qui doit copie le resultat.. (voir)
   return ret;
}
extern Block *currentblock;// def in lg.ypp
void ListOfInst::Add(const C_F0 & ins) {
    if( (!ins.Empty()) ) {
        if( verbosity > 9999 )
            cout << " Add " << n << " " << TheCurrentLine << endl;
        if (n%nx==0){
            Expression   *  l = new Expression [n+nx];
            int * ln =  new int [n+nx];
            int * lsd =  new int [n+nx];
            for (int i=0;i<n;i++) {
                l[i]=list[i];
                ln[i]=linenumber[i];
                lsd[i]=lsldel[i];
            }
            delete [] list;
            delete [] linenumber;
            delete [] lsldel;
            list =l;
            linenumber=ln;
            lsldel=lsd;
        }
        throwassert(list);
        linenumber[n]= TheCurrentLine;
        lsldel[n]=currentblock->nIdWithDelete();
        list[n++] = ins;
    }
}

/// <<ListOfInst::operator()>> Iteratively calls each item in the local array #list of type #Expression

AnyType ListOfInst::operator()(Stack s) const {
    AnyType r;
    int i;
    double s0=CPUtime(),s1=s0,ss0=s0;
    StackOfPtr2Free * sptr = WhereStackOfPtr2Free(s);
    try { // modif FH oct 2006
	for (i=0;i<n;i++)
	{
	    TheCurrentLine=linenumber[i]  ;
            r=(*list[i])(s);
	    sptr->clean(); // modif FH mars 2006  clean Ptr
	    s1=CPUtime();
	    if (showCPU)
             cout << " CPU: "<< i <<" " << linenumber[i] <<  ": " << s1-s0 << "s" << " " << s1-ss0 << "s" << " / " << " " <<lsldel[i] << " " << mpirank << endl;
	    s0=CPUtime();
	}
        if(atclose && atclose->n) {
            if(verbosity>99999 ) cout << " ListOfInst::atclose()  " << n << " " << atclose->n << " // " << lsldel[n-1] << endl;
            atclose->eval(s,atclose->n);}// Add for sep 2016  FH
    }
    catch( E_exception & e)
    {
        if(verbosity>999) cout << " catch E_exception " << i << " " << lsldel[i]  << endl;
        if(atclose) {atclose->eval(s,lsldel[i]);}// Add sep 2016  FH for clean init varaible
	if (e.type() != E_exception::e_return)
	    sptr->clean(); // pour ne pas detruire la valeur retourne  ...  FH  jan 2007
	throw; // rethow
    }
    catch(...)
    {
        if(verbosity>999) cout << " catch ....  " << i << " " << lsldel[i]  << endl;
        if(atclose) {atclose->eval(s,lsldel[i]);}
	sptr->clean();
	throw;
    }
    return r;}

void ShowDebugStack()
 {
   if (mpisize)
   cerr << "  current line = " << TheCurrentLine
        << " mpirank " << mpirank << " / " << mpisize <<endl;
   else
   cerr << "  current line = " << TheCurrentLine  << endl;
  if(debugstack)
      for (int i=0; i<debugstack->size(); ++i)
     {

        cerr << " call " << debugstack->at(i).first->name<< "  at  line "
             <<debugstack->at(i).second << endl;
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
    return Stack_Ptr<AnyType>(s,ParamPtrOffset)[i];
  }
   E_F0para(int ii) : i(ii){}
};

Routine::Routine(aType tf,aType tr,const char * iden,  ListOfId *l,Block * & cb)
    : OneOperator(tr,l),offset(cb->OffSet(sizeof(void*))),
     tfunc(tf),tret(tr),name(iden),param(*l),
      currentblock(new Block(cb)),ins(0)//,clean(0)
     {
       delete l;  // add  FH 24032005 (trap )
       cb = currentblock;
       for (size_t i=0;i<param.size();i++)
       {
           currentblock->NewID(param[i].r,param[i].id,C_F0(new E_F0para(i),// modif FH 2007
							   param[i].r),
							   !param[i].ref);
       }
     }
   Block * Routine::Set(CListOfInst instrs)
       {
           instrs.setclose(Block::snewclose(currentblock));
           ins=instrs;
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
       if(!map_type[l[i].type->name()] )
	 {
	     if(mpirank==0)
	       {
		   cerr << " missing ff type: '" <<l[i].type->name() << "'   "<< map_type.size()  <<  "\n";
		   cerr << "i= " << i << "\n";
	       }
	   InternalError(" missing type ");
	   assert(map_type[l[i].type->name()]);
	 }
       e[i]= map_type[l[i].type->name()]->CastTo(ce);
       k++;
       }
  }

 if (!named_parameter) return;

  if ((size_t) k!=  named_parameter->size())
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
    if ( n && mpirank==0) {
    cerr << " The named parameter can be " << endl;
    for (int i=0;i<n;i++)
       cerr << "\t" << l[i].name << " =  <" << l[i].type->name() << ">\n";
    }
    CompileError("Unused named parameter");
   }
}


//  change FH to bluid .dll

void lgerror (const char* s)
  {
      if(mpirank==0)
	{
	    cerr << endl;
	    cerr <<" Error line number " << zzzfff->lineno() << ", in file " << zzzfff->filename()
	    <<", before  token " <<zzzfff->YYText() << endl
	    << s << endl;
	}
      throw(ErrorCompile(s,zzzfff->lineno(),zzzfff->YYText() ));
  }



 C_F0 ForAll(Block *cb,ListOfId * id,C_F0  m)
{

//Block::open(cb); // new block
     //  decl variable
    if(verbosity>1000)
     cout << "InitAutoLoop ::: " <<  id->size()<< " type=" << *(m.left()) << endl;

     ffassert(id->size()<4);
     aType t=m.left() ;
     ffassert(id->size()<=0 || t->typev);
     ffassert(id->size()<=1 || t->typei);
     ffassert(id->size()<=2 || t->typej);
    // missing this king of code atype<T>->SetTypeLoop(atype<string**>(),atype<K*>())  maybe !!!! FH.
     ffassert(id->size()<4);
     // find the size do data
     aType tt[4];
     int k=0;
     if(t->typei) tt[k++]=t->typei;
     if(t->typej) tt[k++]=t->typej;
     if(t->typev) tt[k++]=t->typev;

     for(int i=0;i<k;++i)
     {
     if(verbosity>1000)
     cout << "     aType = " << i << " left " << *tt[i] << " right="<< *tt[i]->right() <<endl;
     }
     AC_F0 args;
     args=0; // reset
     args+=m;
     for(int j=0,i=id->size(); j<id->size() ; ++j)
     {
         --i;
         C_F0 ci=cb->NewVar<LocalVariable>((*id)[i].id,tt[i]);
         C_F0 cv=Find((*id)[i].id);
         args+=cv;
         const LocalVariable *lv = dynamic_cast<LocalVariable*>((E_F0*) cv);
         if(verbosity>1000)
         cout << " new id " << tt[i] << " "<< (*id)[i].id << " "  << " E="
             <<  (Expression)  args[j] << " "
         <<  (Expression)  ci << " " <<args.size()-1  << " ov: " << (lv ? lv->offset: -1) << " " << *cv.left() << endl;
     }
     Expression loop= new PolymorphicLoop(m,args);
    if(verbosity>1000)
     cout << "a type: " << *atype<PolymorphicLoop*>() << " " << loop << endl;

     return C_F0(loop,atype<PolymorphicLoop*>());
}

 C_F0 ForAll(C_F0  cloop,C_F0  inst)
{
    if(verbosity>1000)
    cout << " type cloop " << *cloop.left() << " " << cloop.LeftValue() << " "  << endl;
    const PolymorphicLoop *loop=  dynamic_cast<const PolymorphicLoop *>(cloop.LeftValue());
    ffassert(loop);
    AC_F0 args;
    args=0;
    args+=loop->t;
    args+=cloop;
    C_F0 instt(inst,atype<NothingType>());
    args+=instt;
    return C_F0(TheOperators,"{}",args);
}

void InitLoop()
{
     Dcl_Type<PolymorphicLoop*>(0);

}
