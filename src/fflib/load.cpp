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
// SUMMARY : Load management
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : Frederic Hecht
// E-MAIL  : frederic.hecht@sorbonne-universite.fr

// TODO: remove this block as soon as autoconf is removed from FreeFEM
#ifndef CMAKE
#include <config.h> // needed for HAVE_DLFCN_H
#endif

#include <iostream>
#include <map>
#include <set>
#include "AFunction.hpp"
#include "environment.hpp"
#include "InitFunct.hpp"
using namespace std;
#include "lex.hpp"
#define LOAD 1

// ALH - 11/3/15 - added NOLOAD macro option because no dynamic load available in javascript

#if defined(__INTEL__) || defined(__MWERKS__) || !defined(HAVE_DLFCN_H) || defined(NOLOAD)
#undef LOAD
#endif

#ifdef LOAD
#include <dlfcn.h>
#elif _WIN32
#include <windows.h>
#endif

#include "ffapi.hpp"

set<string> SetLoadFile;

bool load (string ss) {
  // FFCS - do not allow potentially dangerous commands from remote anonymous clients
    static int count =0;
  if(count++==0) SetLoadFile.insert("msh3");
  if (ffapi::protectedservermode() && (ss == "pipe" || ss == "shell")) {
    cerr << "library " << ss << " not allowed in server environment" << endl;
    CompileError("Error load");
    return 0;
  }

  if (SetLoadFile.find(ss) != SetLoadFile.end()) {
    if ((mpirank == 0) && verbosity)
      cout << " (already loaded: " << ss << ")";
  } else {
    SetLoadFile.insert(ss);
    bool ret = false;
    void *handle = 0;
    const int /*nbprefix=2,*/nbsuffix = 2;
    list<string> prefix(ffenvironment["loadpath"]);
    if (prefix.empty()) {
      prefix.push_back("");
      prefix.push_back("./");
    }

    string suffix[nbsuffix];

    suffix[0] = "";
    suffix[1] = ".so";
#ifdef __APPLE__
    suffix[1] = ".dylib";
#elif _WIN32
    suffix[1] = ".dll";
#endif
    int j;
    for (list<string>::const_iterator i = prefix.begin(); i != prefix.end(); ++i)
      for (j = 0; j < nbsuffix; ++j) {
        string s = *i + ss + suffix[j];

#ifdef LOAD
        handle = dlopen(s.c_str(), RTLD_NOW); // RTLD_GLOBAL RTLD_LAZY
        if (verbosity > 9)
          cout << " test dlopen(" << s << ") = " << handle << endl;

        // FFCS - 20/9/11 - print explanation for load errors
        if (verbosity > 9 && !handle)
          cout << "load error was: " << dlerror() << endl;

        ret = handle != 0;
        if (ret) {
          if (verbosity > 1 && (mpirank == 0))
            cout << " (load: dlopen " << s << " " << handle << ")\n";
          callInitsFunct(); // [[file:InitFunct.cpp::callInitsFunct]]
          return handle;
        }
#elif _WIN32
        {
          HINSTANCE mod = LoadLibrary(s.c_str());
          if (verbosity > 9) cout << " test LoadLibrary(" << s << ") = " << mod << endl;
          if (mod == 0) {
            DWORD merr = GetLastError();
            if (verbosity > 19)
              cerr << "\n try loadLibary: " << s << "\n\t fail: " << merr << endl;
          } else {
            if (verbosity && (mpirank == 0))
              cout << "(load: loadLibary " << s << " = " << handle << ")";
            callInitsFunct(); // [[file:InitFunct.cpp::callInitsFunct]]
            return mod;
          }
        }
#elif STATIC_LINKING
        // <<STATIC_LINKING>> Enable statically linked libraries for [[file:~/fflib/Makefile::STATIC_LINKING]] - ALH
        bool ok = false;

        // <<static_load_msh3>> [[file:~/ff/examples++-load/msh3.cpp::dynamic_loading]]
        if (ss == "msh3-old") {
          // [[file:~/ff/examples++-load/msh3.cpp::msh3_Load_Init]]
          void msh3_Load_Init();
          msh3_Load_Init();
          ok = true;
        }

        // <<static_load_medit>> [[file:~/ff/examples++-load/medit.cpp::dynamic_loading]]
        if (ss == "medit") {
          // [[file:~/ff/examples++-load/medit.cpp::medit_Load_Init]]
          void medit_Load_Init();
          medit_Load_Init();
          ok = true;
        }

        if (ok && verbosity && (mpirank == 0))
          cout << " (static load: " << ss << ")";
        return ok;
#else
        if (mpirank == 0) {
          cout << "--------------------------------------- \n" ;
          cout << "  load: sorry no dlopen on this system " << s << " \n" ;
          cout << "--------------------------------------- \n" ;
        }
        CompileError("not load / dlopen on this system");
        return 0;
#endif
      }

    if (mpirank == 0) {
      cerr << "\nLoad error: " << ss << "\n\t fail: " << endl;
      char *error = 0;
#ifndef _WIN32
#ifdef LOAD
      error = dlerror();
      if (error != NULL) {
        cerr << " dlerror : " << error << endl;
      }
#endif
#endif
      cerr << "list prefix: ";
      for (list<string>::const_iterator i = prefix.begin(); i != prefix.end(); ++i)
        cerr << "'" << *i << "' ";
      cerr << "list suffix: '" << suffix[0] << "' , '" << suffix[1] << "' ";

      cerr << endl;
    }
    throw(ErrorLoad(ss.c_str(), zzzfff->lineno(), zzzfff->YYText()));
    //ErrorLoad("Error load");
  }
  return 0;
}
