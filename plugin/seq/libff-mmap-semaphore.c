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

/* clang-format off */
//ff-c++-LIBRARY-dep:
//ff-c++-cpp-dep:
/* clang-format on */

#include "libff-mmap-semaphore.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef __cplusplus
extern "C" {
#endif

long ff_mmap_sem_verb = 1;

void (*libf_HandleError)(const char *, int err) = 0;
static char *newstringcpy(const char *nm) {
  char *p = (char *)malloc(strlen(nm) + 1);

  if (p) strcpy(p, nm);

  return p;
}

static char *newstringcpyadds(const char *nm) {
  char *p = (char *)malloc(strlen(nm) + 2);

  if (p) {
    if (nm[0] != '/') { /* add / */
      strcpy(p, "/");
      strcat(p, nm);
    } else {
      strcpy(p, nm);
    }
  }

  return p;
}

void ffDoError(const char *msg, int err) {
  if (libf_HandleError) {
    (*libf_HandleError)(msg, err);
  }

  printf(" Error libff-mmap-semaphore: %s Err: %d\n", msg, err);
  exit(1);
}

void ffsem_destroy(ff_Psem p) {
  if (ff_mmap_sem_verb > 9) {
    printf("  ** ffsem_destroy %s unlink: %d\n", p->nm, p->creat);
  }

  int err = 0;
  if (p->creat) {
    err = sem_unlink(p->nm);
  }

  if (err == -1) {
    perror("ff/sem: sem_unlink");
  }

  if (p->sem) {
    err = sem_close(p->sem);
  }

  if (err == -1) {
    perror("ff/sem: sem_close");
  }

  if (p->nm) {
    free((void *)p->nm);
  }

  p->nm = 0;
  p->sem = 0;
}

void ffsem_destroy_(long *p) { ffsem_destroy(*(ff_Psem *)p); }

void ffsem_init0(ff_Psem p) {
  p->sem = 0;
  p->nm = 0;
  p->creat = 0;
}

void ffsem_init0_(long *p) { ffsem_init0(*(ff_Psem *)p); }

ff_Psem ffsem_malloc( ) {
  ff_Psem p = (ff_Psem)malloc(sizeof(struct FF_P_sem));
  if (p == NULL) printf("malloc failed\n");
  p->sem = 0;
  p->nm = 0;
  p->creat = 0;
  return p;
}

void ffsem_del(ff_Psem p) {
  ffsem_destroy(p);
  free(p);
}

void ffsem_del_(long *p) { ffsem_del(*(ff_Psem *)p); }

void ffsem_init(ff_Psem p, const char *nm, int crea) {
  p->creat = crea;
  p->nm = newstringcpyadds(nm);
  if (crea && p->nm) {
    unlink(p->nm);
    p->sem = sem_open(p->nm, O_CREAT, 0660, 0);
  } else {
    p->sem = sem_open(p->nm, 0, 0, 0);
  }

  if (p->sem == SEM_FAILED) {
    p->sem = 0;
    printf(" err sem open %s \n", p->nm);
    perror("sem_open");
    ffsem_destroy(p);
    ffDoError("Error sem_open", 1001);
  }
}

void ffsem_init_(long *pp, const char *nm, int *crea, int lennm) {
  ff_Psem *p = (ff_Psem *)pp;

  *p = ffsem_malloc( );
  ffsem_init(*p, nm, *crea);
}

long ffsem_post(ff_Psem p) {
  int err = sem_post(p->sem);

  if (err == -1) {
    perror("ff/sem: sem_post");
    ffDoError("sem_post", 1002);
  }

  return err;
}

void ffsem_post_(long *p, long *ret) { *ret = ffsem_post(*(ff_Psem *)p); }

long ffsem_wait(ff_Psem p) {
  int err = sem_wait(p->sem);

  if (err == -1) {
    perror("ff/sem: sem_wait");
    ffDoError("sem_post", 1003);
  }

  return err;
}

void ffsem_wait_(long *p, long *ret) { *ret = ffsem_wait(*(ff_Psem *)p); }

long ffsem_trywait(ff_Psem p) {
  int err = sem_trywait(p->sem);

  if (err == -1) {
    perror("ff/sem: sem_trywait");
    ffDoError("sem_post", 1004);
  }

  return err;
}

void ffsem_trywait_(long *p, long *ret) { *ret = ffsem_trywait(*(ff_Psem *)p); }

ff_Pmmap ffmmap_malloc( ) {
  ff_Pmmap p = (ff_Pmmap)malloc(sizeof(struct FF_P_mmap));
  if (p == NULL) printf("malloc failed\n");
  ffmmap_init0(p);
  return p;
}

void ffmmap_del(ff_Pmmap p) {
  ffmmap_destroy(p);
  free(p);
}

void ffmmap_del_(long *p) { ffmmap_del(*(ff_Pmmap *)p); }

void ffmmap_destroy(ff_Pmmap p) {
  if (ff_mmap_sem_verb > 9) {
    printf("  ** ffmmap_destroy %s len: %lu new: %d\n", p->nm, (unsigned long)p->len, p->isnew);
  }

  if (p->map && munmap(p->map, p->len) == -1) {
    printf(" **Error munmap %s %zu\n", p->nm, p->len);
    perror("munmap");
    ffDoError("munmap", 1005);
  }

  if (p->fd > 0) {
    close(p->fd);
  }

  if (p->isnew) {
    unlink(p->nm);
  }

  if (p->nm) {
    free((void *)p->nm);
  }

  p->len = 0;
  p->fd = 0;
  p->nm = 0;
}

void ffmmap_destroy_(long *p) { ffmmap_destroy(*(ff_Pmmap *)p); }

void ffmmap_init0(ff_Pmmap p) {
  p->len = 0;
  p->nm = 0;
  p->fd = 0;
  p->map = 0;
  p->isnew = 0;
}

void ffmmap_init0_(long *p) { ffmmap_init0(*(ff_Pmmap *)p); }

long ffmmap_msync(ff_Pmmap p, long off, long ln) {
  if (ln == 0) {
    ln = p->len - off;
  }

  return msync((char *)p->map + off, ln, MS_SYNC);
}

void ffmmap_msync_(long *p, int *off, int *ln, long *ret) {
  *ret = ffmmap_msync(*(ff_Pmmap *)p, *off, *ln);
}

void ffmmap_init(ff_Pmmap p, const char *nm, long len) {
  void *addr = 0;

  p->len = len;
  p->nm = newstringcpy(nm);    // shm_unlink
  p->map = 0;
  p->fd = open(p->nm, O_RDWR | O_CREAT, (mode_t)0666);
  if (p->fd == -1) {
    printf(" Error opening file mmap  %s  len =  %zu \n", p->nm, p->len);
    perror("open");
    ffmmap_destroy(p);
    ffDoError("opening mmap", 2001);
  }

  off_t size = lseek(p->fd, 0, SEEK_END);    // seek to end of file
  p->isnew = (size == 0);
  printf(" len %ld size %ld \n", len, size);
  if (size < len) {
    if (ftruncate(p->fd, len) == -1) {
      perror("ftruncate");
      printf("Error ftrucated the file %s  len =  %zu \n", p->nm, p->len);
      ffmmap_destroy(p);
      ffDoError("Error ftrucated ", 2002);
    }
  } else {
    p->len = size;
  }

  p->map = mmap(addr, p->len, PROT_READ | PROT_WRITE, MAP_FILE | MAP_SHARED, p->fd, 0);
  if (p->map == MAP_FAILED) {
    p->map = 0;
    printf("Error mmapping the file %s len = %zu\n", p->nm, p->len);
    ffDoError("Error mmapping ", 2003);
  }
}

void ffmmap_init_(long *pp, const char *nm, int *len, int lennm) {
  ff_Pmmap *p = (ff_Pmmap *)pp;

  *p = ffmmap_malloc( );
  ffmmap_init(*p, nm, *len);
}

long ffmmap_read(ff_Pmmap p, void *pt, size_t ln, long off) {
  if (off < 0 || off + ln > p->len) {
    printf("Fatal Error: ffmmap_read ff mmap out of bound len = %zu < %lu + %ld \n", p->len,
           (unsigned long)ln, off);
    ffDoError(" Error out of bound  ", 2004);
  }

  void *pk = (char *)p->map + off;
  memcpy(pt, pk, ln);
  long *pp = (long *)pt;
  if (ff_mmap_sem_verb > 9) {
    printf(" R %ld %ld %lu %p\n", *pp, off, (unsigned long)ln, pk);
  }

  return ln;
}

void ffmmap_read_(long *p, void *pt, int *ln, int *off, long *ret) {
  *ret = ffmmap_read(*(ff_Pmmap *)p, pt, *ln, *off);
  printf("ffmmap_read_ %ld %f %d\n", *(long *)pt, *(double *)pt, *off);
}

long ffmmap_write(ff_Pmmap p, void *pt, size_t ln, long off) {
  if (off < 0 || off + ln > p->len) {
    printf("Fatal Error: ffmmap_write ff mmap out of bound len = %zu < %lu + %ld \n", p->len,
           (unsigned long)ln, off);
    ffDoError(" Error out of bound  ", 2005);
  }

  void *pk = (char *)p->map + off;
  memcpy(pk, pt, ln);
  long *pp = (long *)pk;
  if (ff_mmap_sem_verb > 9) {
    printf(" W %ld %ld %lu %p\n", *pp, off, (unsigned long)ln, pk);
  }

  return ln;
}

void ffmmap_write_(long *p, void *pt, int *ln, int *off, long *ret) {
  printf("ffmmap_write_ %ld %f %d \n", *(long *)pt, *(double *)pt, *off);
  *ret = ffmmap_write(*(ff_Pmmap *)p, pt, *ln, *off);
}

#ifdef __cplusplus
}
#endif
