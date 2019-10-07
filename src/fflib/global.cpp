/// \file

// -*- Mode : c++ -*-
//
// SUMMARY  :
// USAGE    :
// ORG      :
// AUTHOR   : Frederic Hecht
// E-MAIL   : hecht@ann.jussieu.fr
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
#include <iostream>
#include <cstdio>


namespace ffapi {

    //  void init) ();
    // need #include <iostream>
    // need #include <sstream>
    // need using namespace std;
    std::istream * (*cin)();
    std::ostream *(*cout)();
    std::ostream *(*cerr)();

    // <<mingw32_stdout>> Cannot name these functions identically to the original file pointers under MingW32 (compile
    // error). Impacts [[file:InitFunct.hpp::LOADINITIO]]. Changed from stdxxx_ptr() to ffstdxxx() according to the way FF
    // itself was changed.

    FILE *(*ffstdout)();
    FILE *(*ffstderr)();
    FILE *(*ffstdin)();

    /// Initiate graphical pipe output. I need a separate function for this to warn ffcs to check the corresponding ffglut
    /// magic number

    size_t (*fwriteinit)(const void *ptr, size_t size, size_t nmemb,FILE *stream);

    /// Indicates the begining of a new plot to avoid sending socket control data with each plot item.

    void (*newplot)();

    /// Redefinition of standard system calls

    FILE *(*ff_popen)(const char *command, const char *type);
  int (*ff_pclose)(FILE *stream); // [[file:ffapi.cpp::ff_pclose]]
    size_t (*ff_fwrite)(const void *ptr, size_t size, size_t nmemb,FILE *stream);
    int (*ff_fflush)(FILE *stream);
    int (*ff_ferror)(FILE *stream);
    int (*ff_feof)(FILE *stream);

    // Windows file mode
    // -----------------

    /// Changing file mode needs to be disabled when the file is a TCP socket to FFCS. Since the treatment is different in
    /// FF and in FFLANG executables, they have to be stored in a DLL that changes between these two programs.

    void (*wintextmode)(FILE *f);
    void (*winbinmode)(FILE *f);

    // Transfer basic MPI control
    // --------------------------

    void (*mpi_init)(int &argc, char **& argv);
    void (*mpi_finalize)();

    // Permanent server control
    // ------------------------

    /// if true, FF is considered to be accessible from remote anonymous connections and some commands (like shell
    /// commands) are not allowed.

    bool (*protectedservermode)();

}

// TODO: remove this block as soon as autoconf is removed from FreeFem++
#ifndef CMAKE
#include <config.h>
#endif

#include <complex>
#include "AFunction.hpp"
#include "error.hpp"
#include "lex.hpp"
#include "RNM.hpp"
#include <queue>
#include "environment.hpp"
#include "ufunction.hpp"
using namespace std;






#define  FF_GRAPH_PTR_DCL
#include "rgraph.hpp"
#include "fem.hpp"
#include "Mesh3dn.hpp"

#include "HashMatrix.hpp"
#include "SparseLinearSolver.hpp"
#include "MeshPoint.hpp"

 bool  NoGraphicWindow=false;

/// <<verbosity>>
long verbosity = 1;
 long searchMethod = 0; //pichon
long npichon2d=0, npichon3d=0;
long npichon2d1=0, npichon3d1=0;

 FILE *ThePlotStream=0; //  Add for new plot. FH oct 2008

  KN<String> *pkarg;//  for the list of argument  mars 2010
 Map_type_of_map map_type_of_map ; //  to store te type
Map_type_of_map map_pair_of_type ; //  to store te type

 basicForEachType *  typevarreal,  * typevarcomplex;  //  type of real and complex variable

/// <<zzzfff>> see [[file:lex.hpp::mylex]]
mylex *zzzfff;
bool lexdebug;

/// <<plglval>> see [[file:../lglib/lg.ypp::YYSTYPE]] and [[file:../lglib/lg.ypp::yylval]]
#include "lg.tab.hpp"
YYSTYPE *plglval;

 int TheCurrentLine=-1; // unset: by default
//int NbNewVarWithDel =0; // add FH sep 2016 (bof bof global variable not got but hard to set in E_F0 or C_F0
 long mpisize=0,mpirank=0;


   C_F0 *pOne=0,*pZero=0,*pminusOne=0;
// const C_F0 & One(*pOne), &Zero(*pZero);

 Polymorphic * TheOperators=0, //=new Polymorphic(),
             * TheRightOperators=0;//=new Polymorphic();

/// <<Global>> Contains all FreeFem++ language keywords. Declaration in [[file:AFunction.hpp::Global]]

TableOfIdentifier Global;

 long E_Border::Count =0;
 long E_Curve3::Count =0;

/// <<tables_of_identifier>> declared at [[file:AFunction.hpp::tables_of_identifier]]
typedef list<TableOfIdentifier *> ListOfTOfId;
ListOfTOfId tables_of_identifier;

const int AC_F0::MaxSize=1024; // maximal number of parameters



map<const string,basicForEachType *> map_type;
bool showCPU= false;


size_t CodeAlloc::nb=0, CodeAlloc::lg=0,CodeAlloc::nbpx=0,CodeAlloc::chunk=2048;
size_t CodeAlloc::nbt,CodeAlloc::nbdl=0;
CodeAlloc ** CodeAlloc::mem=0;
size_t CodeAlloc::memoryusage=0;
bool CodeAlloc::sort=true;
bool  CodeAlloc::cleanning=false;
bool echo_edp=true; // add F.H of remove script dump

//  add F. Hecht
EnvironmentData  ffenvironment;

basicForEachType *basicForEachType::tnull=0;
E_F0 *E_F0::tnull=0;

long newconvect3=0;// old convect 3d

CodeAlloc *CodeAlloc::tnull=0;

#include <RefCounter.hpp>
RefCounter *RefCounter::tnull=0;
double ff_tgv=1e30;

void InitMeshPoint(void * p)
{
    EF23::MeshPoint*mps=static_cast<EF23::MeshPoint*>(p);
    mps->unset();
}

string *def_solver=0,*def_solver_sym=0, *def_solver_sym_dp=0;
