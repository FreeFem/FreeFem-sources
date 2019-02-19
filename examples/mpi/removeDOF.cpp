// ORIG-DATE: 02/2009
// -*- Mode : c++ -*-
//
// SUMMARY  :  
// USAGE    : LGPL      
// ORG      : LJLL Universite Pierre et Marie Curie, Paris,  FRANCE 
// AUTHOR   : Jacques Morice
// E-MAIL   : jacques.morice@ann.jussieu.fr
//
//ff-c++-LIBRARY-dep:mumps ptscotch scotch metis parmetis  scalapack blas  mpifc  fc mpi  pthread 
//ff-c++-cpp-dep: 

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

 Thank to the ARN ()  FF2A3 grant
 ref:ANR-07-CIS7-002-01 
 */


#include  "ff++.hpp"
static void Load_Init()
{ 
	cout << " load(\"removeDOF\")   is obsoled now,  the function are now in kernel "<< endl; 
	cout <<   "removeDOF ->renumbering\n";
	cout <<   "symmetrizeCSR(A) -> set(A ,sym=1 )\n";
    Global.Add("removeDOF", "(", new removeDOF<std::complex<double>>);
    Global.Add("symmetrizeCSR", "(", new OneOperator1_<long, Matrice_Creuse<double>* >(symmetrizeCSR<double>));

	ExecError("Plugins removeDOF is remove"); 
	exit(0);
}

LOADFUNC(Load_Init);
