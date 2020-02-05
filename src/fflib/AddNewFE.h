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
// to add Finite Element to ff++
class EConstantTypeOfFE :public E_F0
{ 
//  using namespace   Fem2D;
  Fem2D::TypeOfFE * v;
  size_t N;
  bool isconst;
public:
  AnyType operator()(Stack ) const { /*cout << " ()" << v << endl*/;return SetAny<Fem2D::TypeOfFE*>(v);}
    EConstantTypeOfFE( Fem2D::TypeOfFE * o,bool ic=true):v(o),N(v->N),isconst(ic) {assert(v); /*cout << "New constant " << o << endl;*/}
  size_t nbitem() const { return N ;} 
  EConstantTypeOfFE & operator=(Fem2D::TypeOfFE * vv) {
      ffassert( !isconst  && vv && (vv->N == N)); //  same type 
      v = vv;       
      return  *this;
  }
   operator aType () const { assert(v);return atype<Fem2D::TypeOfFE*>();} 
};


struct AddNewFE {  
  AddNewFE (const char * FEname,Fem2D::TypeOfFE* tfe) 
  {
    ffassert(tfe); // check 
    Global.New(FEname, Type_Expr(atype<Fem2D::TypeOfFE*>() ,new  EConstantTypeOfFE(tfe)));
  }
};

// 3d volume case

class EConstantTypeOfFE3 :public E_F0
{ public:
    //  using namespace   Fem2D;
    typedef Fem2D::TypeOfFE3 * T;
    T  v;
public:
    AnyType operator()(Stack ) const { /*cout << " ()" << v << endl*/;return SetAny<T>(v);}
    EConstantTypeOfFE3( T o):v(o) { /*cout << "New constant " << o << endl;*/}
    size_t nbitem() const { assert(v);
    if(verbosity > 2)
        cout << " nb item = " << v->N << endl;
    return v->N ;} 
    operator aType () const { return atype<T>();} 
};

Type_Expr CConstantTFE3(const EConstantTypeOfFE3::T & v);


// 3d surface case

class EConstantTypeOfFES :public E_F0
{ public:
    //  using namespace   Fem2D;
    typedef Fem2D::TypeOfFES * T;
    T  v;
public:
    AnyType operator()(Stack ) const { /*cout << " ()" << v << endl*/;return SetAny<T>(v);}
    EConstantTypeOfFES( T o):v(o) { /*cout << "New constant " << o << endl;*/}
    size_t nbitem() const { assert(v);
        if(verbosity > 2)
            cout << " nb item = " << v->N << endl;
        return v->N ;}
    operator aType () const { return atype<T>();}
};


Type_Expr CConstantTFES(const EConstantTypeOfFES::T & v);


// 3d curve case

class EConstantTypeOfFEL :public E_F0
{ public:
    //  using namespace   Fem2D;
    typedef Fem2D::TypeOfFEL * T;
    T  v;
public:
    AnyType operator()(Stack ) const { /*cout << " ()" << v << endl*/;return SetAny<T>(v);}
    EConstantTypeOfFEL( T o):v(o) { /*cout << "New constant " << o << endl;*/}
    size_t nbitem() const { assert(v);
        if(verbosity > 2)
            cout << " nb item = " << v->N << endl;
        return v->N ;}
    operator aType () const { return atype<T>();}
};


Type_Expr CConstantTFEL(const EConstantTypeOfFEL::T & v);

/*
class EConstantTypeOfFE3 :public E_F0
{ 
    //  using namespace   Fem2D;
    Fem2D::TypeOfFE3 * v;
    size_t N;
    bool isconst;
public:
    AnyType operator()(Stack ) const { return SetAny<Fem2D::TypeOfFE3*>(v);}
    EConstantTypeOfFE3( Fem2D::TypeOfFE3 * o,bool ic=true):v(o),N(v->N),isconst(ic) {assert(v); //cout << "New constant " << o << endl;}
    size_t nbitem() const { return N ;} 
    EConstantTypeOfFE3 & operator=(Fem2D::TypeOfFE3 * vv) {
	ffassert( !isconst  && vv && (vv->N == N)); //  same type 
	v = vv;       
    }
    operator aType () const { assert(v);return atype<Fem2D::TypeOfFE3*>();} 
};
*/
extern map<TypeOfFE *,TypeOfFE3 *> TEF2dto3d;
TypeOfFE * FindFE2(const char * s);

struct AddNewFE3 {  
    AddNewFE3 (const char * FEname,Fem2D::TypeOfFE3* tfe,const char * FEname2=0) 
    {
      ffassert(tfe); // check 
      Global.New(FEname, Type_Expr(atype<Fem2D::TypeOfFE3*>() ,new  EConstantTypeOfFE3(tfe)));
      if(FEname2 && strlen(FEname2))
	  TEF2dto3d[FindFE2(FEname2)]=tfe;
    }
};
extern map<TypeOfFE *,TypeOfFES *> TEF2dtoS;
struct AddNewFES {
    AddNewFES (const char * FEname,Fem2D::TypeOfFES* tfe,const char * FEname2=0)
    {
        ffassert(tfe); // check
        Global.New(FEname, Type_Expr(atype<Fem2D::TypeOfFES*>() ,new  EConstantTypeOfFES(tfe)));
        if(FEname2 && strlen(FEname2))
        TEF2dtoS[FindFE2(FEname2)]=tfe;
    }
};
