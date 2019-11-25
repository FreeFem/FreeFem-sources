#if defined(AIX)
#define wreadmtc wreadmtc
#define userread userread
#define nod2dom nod2dom
#define roscal roscal
#define coscal coscal
#define csrcsc csrcsc
#define aplb aplb
#define expnddom expnddom
#elif defined(BGL)
#define wreadmtc wreadmtc
#define userread userread
#define nod2dom nod2dom
#define roscal roscal
#define coscal coscal
#define csrcsc csrcsc
#define aplb aplb
#define expnddom expnddom
#elif defined(LINUX)
#define wreadmtc wreadmtc_
#define userread userread_
#define nod2dom nod2dom_
#define roscal roscal_
#define coscal coscal_
#define csrcsc csrcsc_
#define aplb aplb_
#define expnddom expnddom_
#elif defined(CRAY)
#define wreadmtc WREADMTC
#define userread USERREAD
#define nod2dom NOD2DOM
#define roscal ROSCAL
#define coscal COSCAL
#define csrcsc CSRCSC
#define aplb APLB
#define expnddom EXPNDDOM
#else
#define wreadmtc wreadmtc_
#define userread userread_
#define nod2dom nod2dom_
#define roscal roscal_
#define coscal coscal_
#define csrcsc csrcsc_
#define aplb aplb_
#define expnddom expnddom_
#endif

/* for outer and inner functions */
extern void setpar(char *filename, char *matrix, int *iov, int *scale, int *unsym, int *method,
                   PrePar prepar, IterPar ipar, DistMatrix dm);

int assignprecon(char *precon_str, DistMatrix dm);

void set_def_params(PrePar prepar, IterPar ipar);

extern void userread(char *fname, int *len, int *job, int *n, int *nnz, double *a, int *ja, int *ia,
                     int *nrhs, double *rhs, int *ierr);
extern void wreadmtc(int *nmax, int *nzmax, int *job, char *fname, int *len, double *a, int *ja,
                     int *ia, double *rhs, int *nrhs, char *guesol, int *nrow, int *ncol, int *nnz,
                     char *title, char *key, char *type, int *ierr);
extern void nod2dom(int *n, int *ndom, int *nodes, int *pointer, int *map, int *mapsize, int *ier);
extern void roscal(int *n, int *job, int *nrm, double *a, int *ja, int *ia, double *diag, double *b,
                   int *jb, int *ib, int *ierr);
extern void coscal(int *n, int *job, int *nrm, double *a, int *ja, int *ia, double *diag, double *b,
                   int *jb, int *ib, int *ierr);
extern void csrcsc(int *n, int *job, int *ipos, double *a, int *ja, int *ia, double *ao, int *jao,
                   int *iao);
extern void aplb(int *nrow, int *ncol, int *job, double *a, int *ja, int *ia, double *b, int *jb,
                 int *ib, double *c, int *jc, int *ic, int *nnzmax, int *iw, int *ierr);
extern void expnddom(int *, int *, int *, int *, int *, int *, int *, int *, int *, int *);

/*-----------------------------------------------------------------------*/
