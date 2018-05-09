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

int P3_Lagrange::nn[10][3] = {
	{0, 0, 0},
	{1, 1, 1},
	{2, 2, 2},
	{1, 1, 2},
	{1, 2, 2},
	{0, 2, 2},
	{0, 0, 2},
	{0, 0, 1},
	{0, 1, 1},
	{0, 1, 2}}
;
int P3_Lagrange::aa[10][3] = {
	{0, 1, 2},
	{0, 1, 2},
	{0, 1, 2},
	{0, 1, 0},
	{0, 0, 1},
	{0, 0, 1},
	{0, 1, 0},
	{0, 1, 0},
	{0, 0, 1},
	{0, 0, 0}}
;
int P3_Lagrange::ff[10] = {6, 6, 6, 2, 2, 2, 2, 2, 2, 1};
int P3_Lagrange::il[10] = {3, 0, 0, 0, 0, 1, 2, 2, 1, 1};
int P3_Lagrange::jl[10] = {0, 3, 0, 2, 1, 0, 0, 1, 2, 1};
int P3_Lagrange::kl[10] = {0, 0, 3, 1, 2, 2, 1, 0, 0, 1};

