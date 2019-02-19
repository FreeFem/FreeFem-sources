/*! \file   elapsed_time.cpp
    \brief  time esurment functions
    \author Atsushi Suzuki, Laboratoire Jacques-Louis Lions
    \date   Jun. 4th 2013
    \date   Jul. 12th 2015
    \date   Nov. 30th 2016
*/

// This file is part of Dissection
// 
// Dissection is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Linking Dissection statically or dynamically with other modules is making
// a combined work based on Disssection. Thus, the terms and conditions of 
// the GNU General Public License cover the whole combination.
//
// As a special exception, the copyright holders of Dissection give you 
// permission to combine Dissection program with free software programs or 
// libraries that are released under the GNU LGPL and with independent modules 
// that communicate with Dissection solely through the Dissection-fortran 
// interface. You may copy and distribute such a system following the terms of 
// the GNU GPL for Dissection and the licenses of the other code concerned, 
// provided that you include the source code of that other code when and as
// the GNU GPL requires distribution of source code and provided that you do 
// not modify the Dissection-fortran interface.
//
// Note that people who make modified versions of Dissection are not obligated 
// to grant this special exception for their modified versions; it is their
// choice whether to do so. The GNU General Public License gives permission to 
// release a modified version without this exception; this exception also makes
// it possible to release a modified version which carries forward this
// exception. If you modify the Dissection-fortran interface, this exception 
// does not apply to your modified version of Dissection, and you must remove 
// this exception when you distribute your modified version.
//
// This exception is an additional permission under section 7 of the GNU 
// General Public License, version 3 ("GPLv3")
//
// Dissection is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Dissection.  If not, see <http://www.gnu.org/licenses/>.
//

#ifndef _elapsed_time_

#ifdef CLOCK_GETTIME

#  ifdef _MSC_VER   // added by Yann Collette
#  include <windows.h>

typedef struct timeval elapsed_t;

#define COPYTIME(a, b) ((a).tv_sec = (b).tv_sec);\
((a).tv_usec = (b).tv_usec)

#  else // ! _MSC_VER == Linux
#include <time.h>
typedef struct timespec elapsed_t;
#define COPYTIME(a, b) ((a).tv_sec = (b).tv_sec);\
((a).tv_nsec = (b).tv_nsec)

#endif

#else   /* #ifdef CLOCK_GETTIME */

#  ifdef GETRUSAGE
#    include <sys/time.h>
#    include <sys/resource.h>
typedef struct rusage elapsed_t;

#define COPYTIME(a, b) ((a).ru_utime.tv_sec = (b).ru_utime.tv_sec); \
((a).ru_utime.tv_usec = (b).ru_utime.tv_usec); \
((a).ru_stime.tv_sec = (b).ru_stime.tv_sec); \
((a).ru_stime.tv_usec = (b).ru_stime.tv_usec) 

#  else /* #ifdef GETTIMEOFDAY */

#ifdef CLOCK // for NEC SX-ACE
#include <time.h>
typedef clock_t elapsed_t;
#define COPYTIME(a, b) (a = b);

#else

#include <sys/time.h>
typedef struct timeval elapsed_t;

#define COPYTIME(a, b) ((a).tv_sec = (b).tv_sec); \
((a).tv_usec = (b).tv_usec)

# endif  /* #ifdef CLOCK */
#  endif /* #ifdef GETRUSAGE */
#endif /* #ifdef CLOCK_GETTIME */

void get_realtime(elapsed_t *tm);
double convert_time(elapsed_t time0, elapsed_t time1);
int convert_sec(elapsed_t t);
int convert_microsec(elapsed_t t);

#define _elapsed_time_
#endif

