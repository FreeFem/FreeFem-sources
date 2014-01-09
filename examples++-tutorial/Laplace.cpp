/*!
 * \file 
 * 
 * \brief C++ version of Laplace.edp test script
 * 
 * 
 * \author Written by Antoine Le Hyaric
 * \author Laboratoire Jacques-Louis Lions
 * \author Université Pierre et Marie Curie-Paris6, UMR 7598, Paris, F-75005 France
 * \author http://www.ljll.math.upmc.fr/lehyaric
 * 
 * \copyright This file is part of Freefem++
 * 
 * \copyright Freefem++ is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2.1 of
 * the License, or (at your option) any later version.
 * 
 * \copyright Freefem++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * \copyright You should have received a copy of the GNU Lesser General Public
 * License along with Freefem++; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
 * 
 * headeralh brief="C++ version of Laplace.edp test script" cpp default=0 dox freefem upmc written
 */

// [[file:~/ff/loc/src/fflib/ffapi.hpp::API]]
#include "../src/fflib/ffapi.hpp"

int main(int argc,char *argv[]){

  // ffapi mention is required for Doxygen to distinguish ffapi::mesh and ::mesh
  ffapi::mesh Th=ffapi::square(10,10);

  fespace Vh(Th,P1);     // P1 FE space
  Vh uh,vh;              // unkown and test function. 
  func f=1;                 //  right hand side function 
  func g=0;                 //  boundary condition function
 
  problem laplace(uh,vh,solver=GMRES,tgv=1e5) =                    //  definion of  the problem 
    int2d(Th)( dx(uh)*dx(vh) + dy(uh)*dy(vh) ) //  bilinear form
    - int2d(Th)( f*vh )                          //  linear form
    + on(1,2,3,4,uh=g) ;                      //  boundary condition form

  laplace; // solve the problem plot(uh); // to see the result
  plot(uh,ps="Laplace.eps",value=true);
}

/*!
 * Local Variables:
 * mode:c++
 * ispell-local-dictionary:"british"
 * coding:utf-8
 * End:
 */
