/// \file
/// ======================================================================
/// Written by Antoine Le Hyaric
/// Laboratoire Jacques-Louis Lions
/// Université Pierre et Marie Curie-Paris6, UMR 7598, Paris, F-75005 France
/// http://www.ljll.math.upmc.fr/lehyaric
/// ======================================================================
/// This file is part of Freefem++
/// 
/// Freefem++ is free software; you can redistribute it and/or modify
/// it under the terms of the GNU Lesser General Public License as
/// published by the Free Software Foundation; either version 2.1 of
/// the License, or (at your option) any later version.
/// 
/// Freefem++  is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU Lesser General Public License for more details.
/// 
/// You should have received a copy of the GNU Lesser General Public
/// License along with Freefem++; if not, write to the Free Software
/// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
/// 02110-1301 USA
/// ======================================================================
/// headeralh cpp freefem start=21/01/10 upmc

#ifndef FFAPI_HPP
#define FFAPI_HPP

// Proposed FreeFem++ Application Programming Interface
// ----------------------------------------------------

// headerfilter

#include <iostream>
#include <sstream>
using namespace std;

namespace ffapi{

  // Redirecting the FF data stream
  // ------------------------------

  // Getting a pointer to FF stdin and stdout enables extra DLLs to
  // use standard IO even when they are redirected (eg through FFCS).

  std::istream *cin();
  std::ostream *cout();
  std::ostream *cerr();

  /// Initiate graphical pipe output. I need a separate function for
  /// this to warn ffcs to check the corresponding ffglut magic number
  size_t fwriteinit(const void *ptr, size_t size, size_t nmemb,FILE *stream);

  /// Indicates the begining of a new plot to avoid sending socket
  /// control data with each plot item.
  void newplot();

  /// Redefinition of standard system calls
  FILE *ff_popen(const char *command, const char *type);
  int ff_pclose(FILE *stream);
  size_t ff_fwrite(const void *ptr, size_t size, size_t nmemb,FILE *stream);
  int ff_fflush(FILE *stream);
  int ff_ferror(FILE *stream);
  int ff_feof(FILE *stream);

  // Windows file mode
  // -----------------

  /// Changing file mode needs to be disabled when the file is a TCP
  /// socket to FFCS. Since the treatment is different in FF and in
  /// FFS executables, they have to be stored in a DLL that changes
  /// between these two programs.
  void wintextmode(FILE *f);
  void winbinmode(FILE *f);

  // Transfer basic MPI control
  // --------------------------

  void mpi_init(int &argc, char **& argv);
  void mpi_finalize();
}

#endif //FFAPI_HPP

/// Local Variables:
/// mode:c++
/// ispell-local-dictionary:"british"
/// coding:utf-8
/// End:
