// SUMMARY  : Threads for Linux and Microsoft Win32
// USAGE    :        
// ORG      : 
// AUTHOR   : Antoine Le Hyaric
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

#include <cassert>
#include <string>
#include <cerrno>
#include <cstring>
#include <iostream>
using namespace std;


#include "ffthreads.hpp"

#ifdef __MINGW32__
#include <process.h>
#else
#include <signal.h>
#endif


Thread::Id Thread::Start(THREADFUNC(f,),Thread::Parm p){
  Id tid;
#ifdef __MINGW32__
  unsigned ThreadId;
  tid = (HANDLE) _beginthreadex(NULL,0,f,p,0,&ThreadId);
  if(tid == NULL) throw StringErr("Thread::Start: Thread could not be created");
#else
  int R = pthread_create(&tid,NULL,f,p);
  if(R != 0) throw StringErr("Thread::Start: Thread could not be created");
#endif
  return tid;
}

void Thread::Wait(Thread::Id tid){
#ifdef __MINGW32__
  DWORD R = WaitForSingleObject(tid,INFINITE);
  if(R == WAIT_FAILED) throw StringErr("Thread::Wait" " -- Wait failed");
  CloseHandle(tid);
#else
  int R=pthread_join(tid,NULL);
  if(R!=0) throw StringErr("Thread::Wait: Wait failed");
#endif
}

void Thread::Exit(){
#ifdef __MINGW32__
    _endthreadex(0);
#else
    pthread_exit(NULL); // No test: returns void.
#endif
}

void Thread::Kill(Thread::Id tid){
#ifdef __MINGW32__
    if(TerminateThread(tid,0) == 0)
      throw StringErr("Thread::Kill: Thread not killed");
#else
    if(pthread_kill(tid,SIGINT)!=0)
      throw StringErr("Thread::Kill: Thread not killed");
#endif
}

Thread::Id Thread::Current(){
#ifdef __MINGW32__
  return GetCurrentThread();
#else
  return pthread_self();
#endif
}
