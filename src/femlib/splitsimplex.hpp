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

template<class Rd> void SplitSimplex (int N, int &nv, Rd * &P, int &nk, int * &K);
// Add J. Morice for function trunc.
void SplitSurfaceSimplex (int N, int &ntri2, int * &tri);

