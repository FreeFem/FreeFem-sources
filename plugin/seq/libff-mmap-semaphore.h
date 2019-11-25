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
// SUMMARY : ...
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : ...
// E-MAIL  : ...

#ifndef _LIBFF_MMAP_SEMAPHORE_
#define _LIBFF_MMAP_SEMAPHORE_

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>
#include <semaphore.h>

#ifdef __cplusplus
extern "C" {
#endif
extern long ff_mmap_sem_verb;
extern void (*libf_HandleError)(const char *, int err);
struct FF_P_sem {
  sem_t *sem;
  const char *nm;
  int creat;
};
struct FF_P_mmap {
  size_t len;
  const char *nm;
  int fd;
  void *map;
  int isnew;
};
typedef struct FF_P_mmap *ff_Pmmap;
typedef struct FF_P_sem *ff_Psem;
ff_Psem ffsem_malloc( );

void ffsem_del(ff_Psem p);
void ffsem_del_(long *p);
void ffsem_destroy(ff_Psem p);
void ffsem_destroy_(long *p);
void ffsem_init0(ff_Psem p);
void ffsem_init0_(long *p);
void ffsem_init(ff_Psem p, const char *nmm, int crea);
void ffsem_init_(long *p, const char *nm, int *crea, int lennm);
long ffsem_post(ff_Psem p);
void ffsem_post_(long *p, long *ret);
long ffsem_wait(ff_Psem p);
void ffsem_wait_(long *p, long *ret);
long ffsem_trywait(ff_Psem p);
void ffsem_trywait_(long *p, long *ret);

ff_Pmmap ffmmap_malloc( );

void ffmmap_del(ff_Pmmap p);
void ffmmap_del_(long *p);
void ffmmap_destroy(ff_Pmmap p);
void ffmmap_destroy_(long *p);
void ffmmap_init0(ff_Pmmap p);
void ffmmap_init0_(long *p);
long ffmmap_msync(ff_Pmmap p, long off, long ln);
void ffmmap_msync_(long *p, int *off, int *ln, long *ret);
void ffmmap_init(ff_Pmmap p, const char *nmm, long len);
void ffmmap_init_(long *p, const char *nm, int *len, int lennm);
long ffmmap_read(ff_Pmmap p, void *t, size_t n, long off);
void ffmmap_read_(long *p, void *pt, int *ln, int *off, long *ret);
long ffmmap_write(ff_Pmmap p, void *t, size_t n, long off);
void ffmmap_write_(long *p, void *pt, int *ln, int *off, long *ret);

#ifdef __cplusplus
}
#endif

#endif
