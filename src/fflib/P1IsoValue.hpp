/*
 * This file is part of FreeFem++.
 *
 * FreeFem++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * FreeFem++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
 */

//
// P1IsoValue.h
// ff
//
// Created by Frédéric Hecht on 07/03/2014.
//
//

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
int IsoLineK (double *f, Fem2D::R3 *Q, double eps);
int IsoLineK (double *f, Fem2D::R2 *Q, double eps);
int UnderIso (double *f, Fem2D::R3 Q[3][4], double vol[3], double eps);
int UnderIso (double *f, Fem2D::R2 Q[2][3], double area[2], double eps);

#endif	/* defined(__ff__P1IsoValue__) */

