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
    int (*ff_pclose)(FILE *stream);
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

#include <cstdio>
#include <complex>
#include <queue>
#include <error.hpp>
#include "environment.hpp"

#define  FF_GRAPH_PTR_DCL
#include "rgraph.hpp"

void ShowDebugStack(){}

 long verbosity = 1;
 long searchMethod=0; // = 9999; //pichon //PROBABLY BUG : can't compile without it
 bool lockOrientation=true;
 FILE *ThePlotStream=0; //  Add for new plot. FH oct 2008


 int TheCurrentLine=-1; // unset: by default
 long mpisize=0,mpirank=0;

bool showCPU= false;
long npichon2d=0, npichon3d=0;
long npichon2d1=0, npichon3d1=0;

//  add F. Hecht
EnvironmentData  ffenvironment;
// April 2019
using std::ostream;
#include <RefCounter.hpp>
RefCounter *RefCounter::tnull=0;
