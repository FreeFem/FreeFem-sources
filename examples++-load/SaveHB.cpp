//ff-c++-LIBRARY-dep: cxx11
//ff-c++-cpp-dep:

// -*- C++ -*-
// Time-stamp: "2018-07-05 14:34:38 fujiwara"
// $Id$
// Author: fujiwara@acs.i.kyoto-u.ac.jp
//----------------------------------------------------------------------
//  Get Linear Equation from FreeFem++ in the HB (Harwell-Boeing) format
//  Japan SIAM, Nagoya, Sept. 2018
//----------------------------------------------------------------------
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <complex>
#include <type_traits>

#include "ff++.hpp"
#include "AFunction_ext.hpp"

void output_matrix_entries( ofstream &fout, const int N,
			    const double *const ccs_val,
			    const int nnzero, const double *const b )
{
    for(int i = 0; i < nnzero; i++){

	fout << setw(20) << setprecision(12) << setiosflags(ios::scientific)
	     << ccs_val[i];

	if((i+1)%4 == 0)
	    fout << std::endl;
    }
    if(nnzero%4 != 0)
	fout << std::endl;

    for(int i = 0; i < N; i++){

	fout << setw(20) << setprecision(12) << setiosflags(ios::scientific)
	     << b[i];

	if((i+1)%4 == 0)
	    fout << std::endl;
    }
    if(N%4 != 0)
	fout << std::endl;

    return;
}

void output_matrix_entries( ofstream &fout, const int N,
			    const complex<double> *const ccs_val,
			    const int nnzero, const complex<double> *const b )
{
    for(int i = 0; i < nnzero; i++){

	fout << setw(20) << setprecision(12) << setiosflags(ios::scientific)
	     << ccs_val[i].real()
	     << setw(20) << setprecision(12) << setiosflags(ios::scientific)
	     << ccs_val[i].imag();

	if((i+1)%2 == 0)
	    fout << std::endl;
    }
    if(nnzero%2 != 0)
	fout << std::endl;

    for(int i = 0; i < N; i++){

	fout << setw(20) << setprecision(12) << setiosflags(ios::scientific)
	     << b[i].real()
	     << setw(20) << setprecision(12) << setiosflags(ios::scientific)
	     << b[i].imag();

	if((i+1)%2 == 0)
	    fout << std::endl;
    }
    if(N%2 != 0)
	fout << std::endl;

    return;
}

template <typename R>
long SaveHB(std::string *const &hb_filename,
	    Matrice_Creuse<R>  *const &sparse_mat,
	    KN_<R> const &b,
	    std::string *const &hb_title)
{
    MatriceMorse<R> *A = sparse_mat->A->toMatriceMorse();

    const bool isDouble = std::is_floating_point<R>::value;

    if(!A)
	return 1; // error

    // A->M = A->N : size of problem
    // A->nbcoef   : # of non-zero entries in the matrix A
    // A->a[i]     : CRS, entry, 0 <= i < A->nbcoef
    // A->cl[i]    : CRS, column indexes, 0 <= i < A->nbcoef
    // A->lg[r]    : CRS, row starting location, 0 <= r < A->m
    // A->lg[A->M] : # of non-zero elements ( == A->nbcoef)

    // b->N()      : size
    // b[i]        : i-th element

    const int N = A->N;
    const int M = A->M;
    ffassert( N == M );

    std::cout << "SaveHB : # of unknowns = " << N << std::endl;

    const int nnzero = A->lg[N];
    ffassert( nnzero == A->nbcoef );

    std::cout << "SaveHB : # of non-zero entries in A = " << nnzero << std::endl;

    //--------------------------------------------------
    // Transpose from CRS to CCS
    //--------------------------------------------------
    int *col_ptr = new int [N+1];
    for(int i = 0; i <= N; i++)    	// initialize
	col_ptr[i] = 0;
    for(int i = 0; i < nnzero; i++)    	// count
	col_ptr[ A->cl[i]+1 ] += 1;
    for(int i = 1; i <= N; i++)		// accumulate
	col_ptr[i] += col_ptr[i-1];

    int *row_ind = new int [ col_ptr[N] ];
    R *ccs_val = new R [ col_ptr[N] ];

    int *count = new int [N];
    for(int i = 0; i < N; i++)
	count[i] = 0;

    for(int i = 0; i < N; i++)
	for(int ij = A->lg[i]; ij < A->lg[i+1]; ij++){

	    int c = A->cl[ij];

	    row_ind[ col_ptr[c] + count[c] ] = i;
	    ccs_val[ col_ptr[c] + count[c] ] = A->a[ij];
	    count[c] += 1;
	}

    delete [] count;

    // from C index manner to FORTRAN index manner
    for(int i = 0; i <= N; i++)
	col_ptr[i]++;
    for(int i = 0; i < nnzero; i++)
	row_ind[i]++;

    //--------------------------------------------------
    // File
    //--------------------------------------------------
    std::ofstream fout(hb_filename->c_str());
    if(!fout){
	std::cerr << "Cannot open the file : "  << hb_filename->c_str() << std::endl;
	std::exit(1);
    }

    std::cout << "SaveHB : filename = " << hb_filename->c_str() << std::endl;

    //--------------------------------------------------
    // Header of the HB format
    //--------------------------------------------------
    const int bufsize = 1024;
    char buf[bufsize];

    // Line 1
    const int HB_TITLE_LENGTH = 72;
    memset(buf, '\0', bufsize);
    strncpy(buf, hb_title->c_str(), HB_TITLE_LENGTH-1);
    std::cout << "SaveHB : title = " << buf << std::endl;

    for(int i = 0; i < HB_TITLE_LENGTH - strlen(hb_title->c_str()); i++)
	strcat(buf, " ");

    strcat(buf, "     KEY");
    fout  << buf << std::endl;

    // Line 2
    int ptrcrd = (N+1)/8;  // #lines for pointers
    if((N+1)%8 != 0)
	ptrcrd += 1;

    int indcrd = nnzero/8;  // #lines for row (or variable) indices
    if(nnzero%8 != 0)
	indcrd += 1;

    int valcrd = (isDouble)? nnzero/4: (2*nnzero)/4; // #lines for numerical values
    if(nnzero%4 != 0)
	valcrd += 1;

    int rhscrd = (isDouble)? N/4: (2*N/4);  // #lines for right-hand sides
    if(N%4 != 0)
	rhscrd += 1;

    // #data lines excluding header
    int totcrd = ptrcrd + indcrd + valcrd + rhscrd;
    sprintf(buf, "%14d%14d%14d%14d%14d", totcrd, ptrcrd, indcrd, valcrd, rhscrd);
    fout << buf << std::endl;

    // Line 3 : MXTYPE, NROW, NCOL, NNZERO, NELTVL
    const int neltvl = 0; // #elemental matrix entries, assembled matrix=0
    if( isDouble )
	sprintf(buf, "RUA           %14d%14d%14d%14d", N, N, nnzero,neltvl);  // double
    else
	sprintf(buf, "CUA           %14d%14d%14d%14d", N, N, nnzero,neltvl);  // complex<double>
    fout << buf << std::endl;

    // Line 4
    //   PTRFMT : Format for Pointers
    //   INDFMT : Format for row (or variable) indices
    //   VALFMT : Format for numerical values of coefficient matrix
    //   RHSFMT : Format for numerical values of right-hand side
    fout <<  "(8I10)          (8I10)          (4E20.12)           (4E20.12)" << std::endl;

    // Line 5
    //  F : Right-hand side type (full storage)
    //  NRHS : #right-hand sides
    //  NRHSIX : #row indices
    fout<< "F                          1             0" << std::endl;

    //--------------------------------------------------
    // Body
    //--------------------------------------------------
    for(int i = 0; i <= N; i++){
	fout << setw(10) << col_ptr[i];
	if((i+1)%8 == 0)
	    fout << std::endl;
    }
    if((N+1)%8 != 0)
	fout << std::endl;

    for(int i = 0; i < nnzero; i++){
	fout << setw(10) << row_ind[i];
	if((i+1)%8 == 0)
	    fout << std::endl;
    }
    if(nnzero%8 != 0)
	fout << std::endl;

    output_matrix_entries( fout, N, ccs_val, nnzero, b );

    fout.close();
    delete [] col_ptr;  delete [] row_ind; delete [] ccs_val;

    delete A;
    return 0;
}
    
static void Load_Init () {    // le constructeur qui ajoute la fonction "splitmesh3"  a freefem++
    // if (verbosity)
    if (verbosity) {cout << " load: SaveHB  " << endl;}
    
	Global.Add("SaveHB", "(", new OneOperator4_<long,string *,Matrice_Creuse<double>*,KN_<double>,string * >(SaveHB));
	Global.Add("SaveHB", "(", new OneOperator4_<long,string *,Matrice_Creuse< std::complex<double> >*,KN_< std::complex<double> >,string * >(SaveHB));
    }
LOADFUNC(Load_Init)
//----------------------------------------------------------------------
// End of file
//----------------------------------------------------------------------
