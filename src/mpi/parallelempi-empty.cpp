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

// empty parallele interface if no MPI to build dll
// with or without mpi..
// #include "parallelepmi.hpp"
extern void (*initparallele)(int &argc, char ** &argv);
extern void (*init_lgparallele)();
extern void (*end_parallele)();

void init_ptr_parallelepmi ();
void init_ptr_parallelepmi () {};

