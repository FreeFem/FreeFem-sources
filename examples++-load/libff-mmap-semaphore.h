
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>
#include <semaphore.h>
extern long ff_mmap_sem_verb;
#ifdef __cplusplus
extern "C" {
#endif

extern void (*libf_HandleError)(const char *,int err);
struct  FF_P_sem {
    sem_t *sem;
    const char *nm;
    int creat;
};

struct FF_P_mmap {
    size_t len;
    const char *nm;
    int fd;
    void * map;
    int isnew;
};
typedef  struct FF_P_mmap *ff_Pmmap ;
typedef  struct FF_P_sem *ff_Psem ;
ff_Psem ffsem_malloc();
void ffsem_del(ff_Psem p);
void ffsem_destroy(ff_Psem p);
void ffsem_init0(ff_Psem p);
void ffsem_init(ff_Psem p,const char *  nmm,int crea);
long ffsem_post(ff_Psem p);
long ffsem_wait(ff_Psem p);
long ffsem_trywait(ff_Psem p);

ff_Pmmap ffmmap_malloc();
void ffmmap_del(ff_Pmmap p);
void ffmmap_destroy(ff_Pmmap p);
void ffmmap_init0(ff_Pmmap p);
long ffmmap_msync(ff_Pmmap p,long off, long ln);
void ffmmap_init(ff_Pmmap p,const char *  nmm,long len);
long ffmmap_read(ff_Pmmap p,void *t,size_t n,size_t off);
long ffmmap_write(ff_Pmmap p,void *t,size_t n,size_t off);

#ifdef __cplusplus
}
#endif
