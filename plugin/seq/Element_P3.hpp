/****************************************************************************/
/* This file is part of FreeFEM.                                            */
/*                                                                          */
/* FreeFEM is free software: you can redistribute it and/or modify          */
/* it under the terms of the GNU Lesser General Public License as           */
/* published by the Free Software Foundation, either version 3 of           */
/* the License, or (at your option) any later version.                      */
/*                                                                          */
/* FreeFEM is distributed in the hope that it will be useful,               */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of           */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            */
/* GNU Lesser General Public License for more details.                      */
/*                                                                          */
/* You should have received a copy of the GNU Lesser General Public License */
/* along with FreeFEM. If not, see <http://www.gnu.org/licenses/>.          */
/****************************************************************************/
// SUMMARY : ...
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : ...
// E-MAIL  : ...

const int TypeOfFE_P3Lagrange::nn[10][3] = {{0, 0, 0}, {1, 1, 1}, {2, 2, 2}, {1, 1, 2}, {1, 2, 2},
                                            {0, 2, 2}, {0, 0, 2}, {0, 0, 1}, {0, 1, 1}, {0, 1, 2}};
const int TypeOfFE_P3Lagrange::aa[10][3] = {{0, 1, 2}, {0, 1, 2}, {0, 1, 2}, {0, 1, 0}, {0, 0, 1},
                                            {0, 0, 1}, {0, 1, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}};
const int TypeOfFE_P3Lagrange::ff[10] = {6, 6, 6, 2, 2, 2, 2, 2, 2, 1};
const int TypeOfFE_P3Lagrange::il[10] = {3, 0, 0, 0, 0, 1, 2, 2, 1, 1};
const int TypeOfFE_P3Lagrange::jl[10] = {0, 3, 0, 2, 1, 0, 0, 1, 2, 1};
const int TypeOfFE_P3Lagrange::kl[10] = {0, 0, 3, 1, 2, 2, 1, 0, 0, 1};
