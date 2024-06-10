/*!
 * \file 
 * 
 * \brief Exposed functions from lgmesh.cpp
 * 
 * 
 * \author Written by Antoine Le Hyaric
 * \author http://www.ljll.math.upmc.fr/lehyaric
 * \author Laboratoire Jacques-Louis Lions
 * \author Université Pierre et Marie Curie-Paris6, UMR 7598, Paris, F-75005 France
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
 * headeralh brief="Exposed functions from lgmesh.cpp" cpp default=0 dox freefem start=21/10/2013 upmc written
 */

  // <<Carre>> [[file:lgmesh.cpp::Carre]]
  Fem2D::Mesh * Carre_(int nx,int ny,Expression fx,Expression fy,Stack stack,int flags,KN_<long> lab,long reg=0,bool rmdup=false);
  Fem2D::Mesh* Carre(int nx,int ny,Expression fx,Expression fy,Stack stack,int flags,KN_<long> lab,long reg=0,bool rmdup=false);

/*!
 * Local Variables:
 * mode:c++
 * ispell-local-dictionary:"british"
 * coding:utf-8
 * End:
 */
