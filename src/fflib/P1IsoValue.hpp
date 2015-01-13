//
//  P1IsoValue.h
//  ff
//
//  Created by Frédéric Hecht on 07/03/2014.
//
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


#ifndef __ff__P1IsoValue__
#define __ff__P1IsoValue__

#include <cmath>
#include <cstdlib>
#include <iostream>


#include "ufunction.hpp"
namespace Fem2D
{
#include "R3.hpp"
}
int IsoLineK(double *f,Fem2D::R3 *Q, double eps);
int IsoLineK(double *f,Fem2D::R2 *Q, double eps);
int UnderIso(double *f,Fem2D::R3 Q[3][4] ,double vol[3],  double eps);
int UnderIso(double *f,Fem2D::R2 Q[2][3] ,double area[2],  double eps);



#endif /* defined(__ff__P1IsoValue__) */
