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
public:
  AnyType operator()(Stack ) const { /*cout << " ()" << v << endl*/;return SetAny<Fem2D::TypeOfFE*>(v);}
  EConstantTypeOfFE( Fem2D::TypeOfFE * o):v(o) { /*cout << "New constant " << o << endl;*/}
  size_t nbitem() const { assert(v);return v->N ;} 
   operator aType () const { return atype<Fem2D::TypeOfFE*>();} 
};


struct AddNewFE {  
  AddNewFE (const char * FEname,Fem2D::TypeOfFE* tfe) 
  {
    ffassert(tfe); // check 
    Global.New(FEname, Type_Expr(atype<Fem2D::TypeOfFE*>() ,new  EConstantTypeOfFE(tfe)));
  }
};
