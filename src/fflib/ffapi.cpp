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

// Proposed FreeFem++ Application Programming Interface
// ----------------------------------------------------

// headerfilter
#include "ffapi.hpp"
#ifndef FFS
#include <cstdio>
#endif
#ifdef FFS
#include "../src/options.hpp"
#endif
#ifndef FFS
#ifdef WIN32
#include <fcntl.h>
#endif
#endif
#ifndef FFS
#ifdef PARALLELE
#include "mpi.h"
#endif
#endif

#include <cstdlib>

// FFCS-specific implementations for the FF API
// --------------------------------------------

/// FFCS defined means that FFCS is being compiled. I am fairly
/// confident that FFCS will not be defined while compiling the
/// original FF.
#ifdef FFS
#include "../src/socket.hpp"
#include "../src/buffer.hpp"
#endif

/// Need to choose a non-zero stream number because FF will check it
/// (as global variable ThePlotStream)
#define FFAPISTREAM 1

/// if FFCS is around, we need to bufferize all communication to avoid
/// mixing up CMD_FFG and CMD_STDOUT messages
#ifdef FFS
void bufferwrite(const char *b,const int l){

  // thank to the buffering, there is only one CMD_FFG tag for multiple
  // visualization data items.
  *serversocket<<CMD_FFG;
  *serversocket<<l;

  // this call contains the socket MAGIC number
  serversocket->bufferedwrite(static_cast<const char*>(b),l);
}

Buffer buffer(NULL,bufferwrite);
#endif

namespace ffapi{

  // Get a pointer to the local cin/cout (which is distinct from
  // ffcs's stdin/stdout under Windows because each DLL owns separate
  // cin/cout objects).
  std::istream *cin(){return &std::cin;}
  std::ostream *cout(){return &std::cout;}
  std::ostream *cerr(){return &std::cerr;}

  void newplot(){
#ifdef FFS
    assert(serversocket);
#endif
  }

  FILE *ff_popen(const char *command, const char *type){
#ifdef FFS
    // this happens right at the begining of FF, so the socket
    // communication must not be started yet (only when actual
    // visualization data needs to be transfered).
    return (FILE*)FFAPISTREAM;
#else
    // need #include <cstdio>
    popen(command,type);
#endif
  }

  int ff_pclose(FILE *stream){
#ifdef FFS
    // nothing to close in FFCS
    return 0;
#else
    pclose(stream);
#endif
  }

  size_t fwriteinit(const void *ptr, size_t size, size_t nmemb,FILE *stream){

    // printf() is useful for debug because it is not redirected through
    // the FFCS socket. But it is asynchronous with cout so it may end up
    // in the middle of the lines checked by test/compare. So deactivate
    // it by default.
#ifdef DEBUG_FFAPI
#ifdef FFS
    printf("debug: ffapi: using TCP sockets\n");
#else
    printf("debug: ffapi: using an anonymous pipe\n");
#endif // FFS
#endif // DEBUG_FFAPI

#ifdef FFS
    // Ask FFCS to analyze the visualization flux header. I could just
    // skip this stage, but it will be useful to check the coherency
    // between FFCS and FF when FF evolves in the future.
    assert(serversocket);
    *serversocket<<CMD_FFGINIT;
#endif
    ff_fwrite(ptr,size,nmemb,stream);
  }

  size_t ff_fwrite(const void *ptr, size_t size, size_t nmemb,FILE *stream){
#ifdef FFS
    // this assert is a way to check that the serversocket pointer
    // exported from the FFCS shared library is a valid one (which has
    // not been always true in the case of Windows DLLs).
    assert(serversocket);

    // we won't make use of the stream, but make sure that the call from
    // FF is coherent with what we know.
    assert(stream==(FILE*)FFAPISTREAM);

    buffer.write(static_cast<const char*>(ptr),size*nmemb);

    // stops the server flux at one precise point (point value expressed
    // during a previous crash while reading server data in the client
    // in visudata.cpp). Use abort() to call the debugger (which can
    // display the call stack and show where the problematic pipe value
    // came from).

    // need #include "../src/options.hpp"
    if(options->AbortFFGDataAt==buffer.getpoint())abort();

#else
    fwrite(ptr,size,nmemb,stream);
#endif
  }

  int ff_fflush(FILE *stream){
#ifdef FFS
    assert(stream==(FILE*)FFAPISTREAM);
    assert(serversocket);

    // we need to flush both the buffer and the socket to avoid a
    // separate callback for flush in the buffer
    buffer.flush();
    serversocket->writeflush();
#else
    fflush(stream);
#endif
  }

  int ff_ferror(FILE *stream){
#ifndef FFS
    return ferror(stream);
#endif
  }

  int ff_feof(FILE *stream){
#ifndef FFS
    return feof(stream);
#endif
  }

  void wintextmode(FILE *f){
#ifndef FFS
#ifdef WIN32
    // need #include <fcntl.h>
    _setmode(fileno(f),O_TEXT);	
#endif
#endif
  }

  void winbinmode(FILE *f){
#ifndef FFS
#ifdef WIN32
    _setmode(fileno(f),O_BINARY);	
#endif
#endif
  }

  void mpi_init(int &argc, char** &argv){
    /// only call MPI_Init() if this has not already been done in
    /// ffcs/src/server.cpp
#ifndef FFS
#ifdef PARALLELE
    // need #include "mpi.h"
    MPI_Init(&argc,&argv);
#endif
#endif
  }

  void mpi_finalize(){
#ifndef FFS
#ifdef PARALLELE
    MPI_Finalize();
#endif
#endif
  }

void init()
{
  ffapi::fwriteinit ;
  ffapi::winbinmode ;
  ffapi::wintextmode ;
  ffapi::mpi_finalize ;
  ffapi::cin ;
  ffapi::cerr ;
  ffapi::cout ;
  ffapi::ff_feof ;
  ffapi::newplot ;
  ffapi::ff_popen ;
  ffapi::mpi_init ;
  ffapi::ff_ferror ;
  ffapi::ff_fflush ;
  ffapi::ff_fwrite ;
  ffapi::ff_pclose ;
;

}
}

/// Local Variables:
/// mode:c++
/// ispell-local-dictionary:"british"
/// coding:utf-8
/// End:
