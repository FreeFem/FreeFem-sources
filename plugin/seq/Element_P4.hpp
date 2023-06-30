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

const int TypeOfFE_P4Lagrange::nn[15][4] = {{0, 0, 0, 0}, {1, 1, 1, 1}, {2, 2, 2, 2}, {1, 1, 1, 2},
                                            {1, 1, 2, 2}, {1, 2, 2, 2}, {0, 2, 2, 2}, {0, 0, 2, 2},
                                            {0, 0, 0, 2}, {0, 0, 0, 1}, {0, 0, 1, 1}, {0, 1, 1, 1},
                                            {0, 1, 2, 2}, {0, 1, 1, 2}, {0, 0, 1, 2}};
const int TypeOfFE_P4Lagrange::aa[15][4] = {{0, 1, 2, 3}, {0, 1, 2, 3}, {0, 1, 2, 3}, {0, 1, 2, 0},
                                            {0, 1, 0, 1}, {0, 0, 1, 2}, {0, 0, 1, 2}, {0, 1, 0, 1},
                                            {0, 1, 2, 0}, {0, 1, 2, 0}, {0, 1, 0, 1}, {0, 0, 1, 2},
                                            {0, 0, 0, 1}, {0, 0, 1, 0}, {0, 1, 0, 0}};
const int TypeOfFE_P4Lagrange::ff[15] = {24, 24, 24, 6, 4, 6, 6, 4, 6, 6, 4, 6, 2, 2, 2};
const int TypeOfFE_P4Lagrange::il[15] = {4, 0, 0, 0, 0, 0, 1, 2, 3, 3, 2, 1, 1, 1, 2};
const int TypeOfFE_P4Lagrange::jl[15] = {0, 4, 0, 3, 2, 1, 0, 0, 0, 1, 2, 3, 1, 2, 1};
const int TypeOfFE_P4Lagrange::kl[15] = {0, 0, 4, 1, 2, 3, 3, 2, 1, 0, 0, 0, 2, 1, 1};

const int TypeOfFE_P4_S::nn[15][4] = {{0, 0, 0, 0}, {1, 1, 1, 1}, {2, 2, 2, 2}, {1, 1, 1, 2},
                                            {1, 1, 2, 2}, {1, 2, 2, 2}, {0, 2, 2, 2}, {0, 0, 2, 2},
                                            {0, 0, 0, 2}, {0, 0, 0, 1}, {0, 0, 1, 1}, {0, 1, 1, 1},
                                            {0, 1, 2, 2}, {0, 1, 1, 2}, {0, 0, 1, 2}};
const int TypeOfFE_P4_S::aa[15][4] = {{0, 1, 2, 3}, {0, 1, 2, 3}, {0, 1, 2, 3}, {0, 1, 2, 0},
                                            {0, 1, 0, 1}, {0, 0, 1, 2}, {0, 0, 1, 2}, {0, 1, 0, 1},
                                            {0, 1, 2, 0}, {0, 1, 2, 0}, {0, 1, 0, 1}, {0, 0, 1, 2},
                                            {0, 0, 0, 1}, {0, 0, 1, 0}, {0, 1, 0, 0}};
const int TypeOfFE_P4_S::ff[15] = {24, 24, 24, 6, 4, 6, 6, 4, 6, 6, 4, 6, 2, 2, 2};
const int TypeOfFE_P4_S::il[15] = {4, 0, 0, 0, 0, 0, 1, 2, 3, 3, 2, 1, 1, 1, 2};
const int TypeOfFE_P4_S::jl[15] = {0, 4, 0, 3, 2, 1, 0, 0, 0, 1, 2, 3, 1, 2, 1};
const int TypeOfFE_P4_S::kl[15] = {0, 0, 4, 1, 2, 3, 3, 2, 1, 0, 0, 0, 2, 1, 1};
