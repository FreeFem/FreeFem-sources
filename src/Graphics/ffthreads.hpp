// SUMMARY  : Threads for Linux and Microsoft Win32
// USAGE    :        
// ORG      : 
// AUTHOR   : Antoine Le Hyaric, Modif F. hecht
// E-MAIL   : lehyaric@ann.jussieu.fr

// This file is part of Freefem++
// 
// Freefem++ is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation; either version 2.1 of the License, or
// (at your option) any later version.
// 
// Freefem++  is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with Freefem++; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

// Antoine Le Hyaric - LJLL Paris 6 - lehyaric@ann.jussieu.fr - 21/10/04

#ifndef THREADS_HPP
#define THREADS_HPP

#include <string>
using namespace std;
#ifdef __MINGW32__
#include <windows.h>
#else
#include <pthread.h>
#endif

// Just check that we are in a known environment (otherwise it may be
// difficult to recognise that the simple cause is an ifdef problem).

class Thread{
public:

#ifdef __MINGW32__
#define THREADFUNC(f,parm) unsigned int (__stdcall f)(Thread::Parm parm)
  typedef LPVOID Parm;
  typedef HANDLE Id;
#else
#define THREADFUNC(f,parm) void* f(Thread::Parm parm)
  typedef void* Parm;
  typedef pthread_t Id;
#endif

  // Mingw is a little puzzled if there are no brackets around
  // __stdcall
  static Id Start(THREADFUNC(f,),Parm p);
  static void Wait(Id tid);
  static void Exit(); // From inside the thread
  static void Kill(Id tid);

  static Id Current();
};
 
struct StringErr {
  const char * s;
  StringErr(const char * ss) : s(ss) {}
  friend ostream & operator<<(ostream & f,const StringErr & se) { return f << se.s << endl;}
};

#endif // THREADS_HPP
