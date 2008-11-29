/*
  The freefem++ arpack interface without arpack ++ 
  aprack++ is do  too much c++ error  with modern  compiler
  but a by tanks  to arpack++ for the idea of this code
   F. Hecht. 
*/
//  Change V[1] -> V[0]
//  Change workev[1] -> workev[0]
#ifndef ARPACKFF_HPP
#define ARPACKFF_HPP
#include "error.hpp"

typedef int integer;
typedef int logical;

#define F77NAME(x) x ## _

#define HIDDEN_HBW ,int,int,int
#define HIDDEN_BW ,int,int
#define HIDDEN_12 ,1,2
#define HIDDEN_112 ,1,1,2


extern "C"
{
  /*
// debug "common" statement.

  struct { 
    integer logfil, ndigit, mgetv0;
    integer msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd;
    integer mnaupd, mnaup2, mnaitr, mneigt, mnapps, mngets, mneupd;
    integer mcaupd, mcaup2, mcaitr, mceigt, mcapps, mcgets, mceupd;
  } F77NAME(debug);

  // add FH pour lib dynamic so win 32 
  struct { 
    integer           nopx, nbx, nrorth, nitref, nrstrt;
    float           tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv,
      tnaupd, tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv,
      tcaupd, tcaup2, tcaitr, tceigh, tcgets, tcapps, tcconv,
      tmvopx, tmvbx, tgetv0, titref, trvec;
    
    
  } F77NAME(timing);

*/
// double precision symmetric routines.

  void F77NAME(dsaupd)(integer *ido, char *bmat, integer *n, char *which,
                       integer *nev, double *tol, double *resid,
                       integer *ncv, double *V, integer *ldv,
                       integer *iparam, integer *ipntr, double *workd,
                       double *workl, integer *lworkl, integer *info HIDDEN_BW);

  void F77NAME(dseupd)(logical *rvec, char *HowMny, logical *select,
                       double *d, double *Z, integer *ldz,
                       double *sigma, char *bmat, integer *n,
                       char *which, integer *nev, double *tol,
                       double *resid, integer *ncv, double *V,
                       integer *ldv, integer *iparam, integer *ipntr,
                       double *workd, double *workl,
                       integer *lworkl, integer *info HIDDEN_HBW );

// double precision nonsymmetric routines.

  void F77NAME(dnaupd)(integer *ido, char *bmat, integer *n, char *which,
                       integer *nev, double *tol, double *resid,
                       integer *ncv, double *V, integer *ldv,
                       integer *iparam, integer *ipntr, double *workd,
                       double *workl, integer *lworkl, integer *info HIDDEN_BW);

  void F77NAME(dneupd)(logical *rvec, char *HowMny, logical *select,
                       double *dr, double *di, double *Z,
                       integer *ldz, double *sigmar,
                       double *sigmai, double *workev,
                       char *bmat, integer *n, char *which,
                       integer *nev, double *tol, double *resid,
                       integer *ncv, double *V, integer *ldv,
                       integer *iparam, integer *ipntr,
                       double *workd, double *workl,
                       integer *lworkl, integer *info HIDDEN_HBW);

// single precision symmetric routines.

  void F77NAME(ssaupd)(integer *ido, char *bmat, integer *n, char *which,
                       integer *nev, float *tol, float *resid,
                       integer *ncv, float *V, integer *ldv,
                       integer *iparam, integer *ipntr, float *workd,
                       float *workl, integer *lworkl, integer *info HIDDEN_BW);

  void F77NAME(sseupd)(logical *rvec, char *HowMny, logical *select,
                       float *d, float *Z, integer *ldz,
                       float *sigma, char *bmat, integer *n,
                       char *which, integer *nev, float *tol,
                       float *resid, integer *ncv, float *V,
                       integer *ldv, integer *iparam, integer *ipntr,
                       float *workd, float *workl,
                       integer *lworkl, integer *info HIDDEN_HBW);

// single precision nonsymmetric routines.

  void F77NAME(snaupd)(integer *ido, char *bmat, integer *n, char *which,
                       integer *nev, float *tol, float *resid,
                       integer *ncv, float *V, integer *ldv,
                       integer *iparam, integer *ipntr, float *workd,
                       float *workl, integer *lworkl, integer *info HIDDEN_BW);

  void F77NAME(sneupd)(logical *rvec, char *HowMny, logical *select,
                       float *dr, float *di, float *Z,
                       integer *ldz, float *sigmar,
                       float *sigmai, float *workev, char *bmat,
                       integer *n, char *which, integer *nev,
                       float *tol, float *resid, integer *ncv,
                       float *V, integer *ldv, integer *iparam,
                       integer *ipntr, float *workd, float *workl,
                       integer *lworkl, integer *info HIDDEN_HBW);


  void F77NAME(cnaupd)(integer *ido, char *bmat, integer *n, char *which,
                       integer *nev, float *tol, complex<float> *resid,
                       integer *ncv, complex<float> *V, integer *ldv,
                       integer *iparam, integer *ipntr, complex<float> *workd,
                       complex<float> *workl, integer *lworkl,
                       float *rwork, integer *info HIDDEN_BW);

  void F77NAME(cneupd)(logical *rvec, char *HowMny, logical *select,
                       complex<float> *d, complex<float> *Z, integer *ldz,
                       complex<float> *sigma, complex<float> *workev,
                       char *bmat, integer *n, char *which, integer *nev,
                       float *tol, complex<float> *resid, integer *ncv,
                       complex<float> *V, integer *ldv, integer *iparam,
                       integer *ipntr, complex<float> *workd,
                       complex<float> *workl, integer *lworkl,
                       float *rwork, integer *info HIDDEN_HBW);

// double precision complex routines.

  void F77NAME(znaupd)(integer *ido, char *bmat, integer *n, char *which,
                       integer *nev, double *tol, complex<double> *resid,
                       integer *ncv, complex<double> *V, integer *ldv,
                       integer *iparam, integer *ipntr, complex<double> *workd,
                       complex<double> *workl, integer *lworkl,
                       double *rwork, integer *info HIDDEN_BW);

  void F77NAME(zneupd)(logical *rvec, char *HowMny, logical *select,
                       complex<double> *d, complex<double> *Z, integer *ldz,
                       complex<double> *sigma, complex<double> *workev,
                       char *bmat, integer *n, char *which, integer *nev,
                       double *tol, complex<double> *resid, integer *ncv,
                       complex<double> *V, integer *ldv, integer *iparam,
                       integer *ipntr, complex<double> *workd,
                       complex<double> *workl, integer *lworkl,
                       double *rwork, integer *info HIDDEN_HBW);

}

inline void saupp(int& ido, char bmat, int n, char* which, int nev,
                  double& tol, double resid[], int ncv, double V[],
                  int ldv, int iparam[], int ipntr[], double workd[],
                  double workl[], int lworkl, int& info)

{

  F77NAME(dsaupd)(&ido, &bmat, &n, which, &nev, &tol, resid, &ncv,
                  &V[0], &ldv, &iparam[1], &ipntr[1], &workd[1], &workl[1],
                  &lworkl, &info HIDDEN_12);

} // saupp (double).

inline void saupp(int& ido, char bmat, int n, char* which, int nev,
                  float& tol, float resid[], int ncv, float V[],
                  int ldv, int iparam[], int ipntr[], float workd[],
                  float workl[], int lworkl, int& info)

{

  F77NAME(ssaupd)(&ido, &bmat, &n, which, &nev, &tol, resid, &ncv,
                  &V[0], &ldv, &iparam[1], &ipntr[1], &workd[1], &workl[1],
                  &lworkl, &info HIDDEN_12);

} // saupp (float).


inline void seupp(bool rvec, char HowMny, double d[], double Z[],
                  int ldz, double sigma, char bmat, int n,
                  char* which, int nev, double tol, double resid[],
                  int ncv, double V[], int ldv, int iparam[],
                  int ipntr[], double workd[], double workl[],
                  int lworkl, int& info)

{

  int      irvec;
  logical* iselect;
  double*  iZ;

  irvec   = (int) rvec;
  iselect = new logical[ncv];
  iZ = (Z == NULL) ? &V[0] : Z;

  F77NAME(dseupd)(&irvec, &HowMny, iselect, d, iZ, &ldz, &sigma, &bmat,
                  &n, which, &nev, &tol, resid, &ncv, &V[0], &ldv, &iparam[1],
                  &ipntr[1], &workd[1], &workl[1], &lworkl, &info HIDDEN_112);

  delete[] iselect;

} // seupp (double).
inline void seupp(bool rvec, char HowMny, float d[], float Z[],
                  int ldz, float sigma, char bmat, int n,
                  char* which, int nev, float tol, float resid[],
                  int ncv, float V[], int ldv, int iparam[],
                  int ipntr[], float workd[], float workl[],
                  int lworkl, int& info)
{

  int      irvec;
  logical* iselect;
  float*   iZ;

  irvec   = (int) rvec;
  iselect = new logical[ncv];
  iZ = (Z == NULL) ? &V[0] : Z;

  F77NAME(sseupd)(&irvec, &HowMny, iselect, d, iZ, &ldz, &sigma, &bmat,
                  &n, which, &nev, &tol, resid, &ncv, &V[0], &ldv, &iparam[1],
                  &ipntr[1], &workd[1], &workl[1], &lworkl, &info HIDDEN_112 );

  delete[] iselect;

} // seupp (float).

inline void naupp(int& ido, char bmat, int n, char* which, int nev,
                  double& tol, double resid[], int ncv, double V[],
                  int ldv, int iparam[], int ipntr[], double workd[],
                  double workl[], int lworkl, int& info)
{

  F77NAME(dnaupd)(&ido, &bmat, &n, which, &nev, &tol, resid, &ncv,
                  &V[0], &ldv, &iparam[1], &ipntr[1], &workd[1], &workl[1],
                  &lworkl, &info HIDDEN_12);

} // naupp (double).

inline void naupp(int& ido, char bmat, int n, char* which, int nev,
                  float& tol, float resid[], int ncv, float V[],
                  int ldv, int iparam[], int ipntr[], float workd[],
                  float workl[], int lworkl, int& info)

{

  F77NAME(snaupd)(&ido, &bmat, &n, which, &nev, &tol, resid, &ncv,
                  &V[0], &ldv, &iparam[1], &ipntr[1], &workd[1], &workl[1],
                  &lworkl, &info HIDDEN_12 );

} // naupp (float).



inline void caupp(int& ido, char bmat, int n, char* which, int nev,
                  double& tol, complex<double> resid[], int ncv,
                  complex<double> V[], int ldv, int iparam[], int ipntr[],
                  complex<double> workd[], complex<double> workl[],
                  int lworkl, double rwork[], int& info)

{

  F77NAME(znaupd)(&ido, &bmat, &n, which, &nev, &tol, resid, &ncv,
                  &V[0], &ldv, &iparam[1], &ipntr[1], &workd[1],
                  &workl[1], &lworkl, &rwork[1], &info HIDDEN_12 );

} // caupp 

inline void caupp(int& ido, char bmat, int n, char* which, int nev,
                  float& tol, complex<float> resid[], int ncv,
                  complex<float> V[], int ldv, int iparam[], int ipntr[],
                  complex<float> workd[], complex<float> workl[],
                  int lworkl, float rwork[], int& info)
{

  F77NAME(cnaupd)(&ido, &bmat, &n, which, &nev, &tol, resid, &ncv,
                  &V[0], &ldv, &iparam[1], &ipntr[1], &workd[1],
                  &workl[1], &lworkl, &rwork[1], &info HIDDEN_12);

} // caupp 

inline void ceupp(bool rvec, char HowMny, complex<double> d[],
		 complex<double> Z[], int ldz, complex<double> sigma,
		 complex<double> workev[], char bmat, int n, char* which,
		 int nev, double tol, complex<double> resid[], int ncv,
		 complex<double> V[], int ldv, int iparam[], int ipntr[],
		 complex<double> workd[], complex<double> workl[],
		 int lworkl, double rwork[], int& info)
{

  int                irvec;
  logical*           iselect;
  complex<double>* iZ;

  irvec   = (int) rvec;
  iselect = new logical[ncv];
  iZ = (Z == NULL) ? &V[0] : Z;

  F77NAME(zneupd)(&irvec, &HowMny, iselect, d, iZ, &ldz, &sigma,
                  &workev[0], &bmat, &n, which, &nev, &tol, resid,
                  &ncv, &V[0], &ldv, &iparam[1], &ipntr[1],
                  &workd[1], &workl[1], &lworkl, &rwork[1], &info HIDDEN_112);

  delete[] iselect;

} // ceupp (complex<double>).

inline void ceupp(bool rvec, char HowMny, complex<float> d[],
                  complex<float> Z[], int ldz, complex<float> sigma,
                  complex<float> workev[], char bmat, int n, char* which,
                  int nev, float tol, complex<float> resid[], int ncv,
                  complex<float> V[], int ldv, int iparam[], int ipntr[],
                  complex<float> workd[], complex<float> workl[],
                  int lworkl, float rwork[], int& info)

{

  int               irvec;
  logical*          iselect;
  complex<float>* iZ;

  irvec   = (int) rvec;
  iselect = new logical[ncv];
  iZ = (Z == NULL) ? &V[0] : Z;

  F77NAME(cneupd)(&irvec, &HowMny, iselect, d, iZ, &ldz, &sigma,
                  &workev[0], &bmat, &n, which, &nev, &tol, resid,
                  &ncv, &V[0], &ldv, &iparam[1], &ipntr[1],
                  &workd[1], &workl[1], &lworkl, &rwork[1], &info HIDDEN_112);

  delete[] iselect;

} // ceupp (complex<float>).

inline void neupp(bool rvec, char HowMny, double dr[],
		 double di[], double Z[], int ldz, double sigmar,
		 double sigmai, double workv[], char bmat, int n,
		 char* which, int nev, double tol, double resid[],
		 int ncv, double V[], int ldv, int iparam[],
		 int ipntr[], double workd[], double workl[],
		 int lworkl, int& info)
{

  int      irvec;
  logical* iselect;
  double*  iZ;

  irvec   = (int) rvec;
  iselect = new logical[ncv];
  iZ = (Z == NULL) ? &V[0] : Z;

  F77NAME(dneupd)(&irvec, &HowMny, iselect, dr, di, iZ, &ldz, &sigmar,
                  &sigmai, &workv[1], &bmat, &n, which, &nev, &tol,
                  resid, &ncv, &V[0], &ldv, &iparam[1], &ipntr[1],
                  &workd[1], &workl[1], &lworkl, &info HIDDEN_112);

  delete[] iselect;

} // neupp (double).                                                                                  

inline void neupp(bool rvec, char HowMny, float dr[],
                  float di[], float Z[], int ldz, float sigmar,
                  float sigmai, float workv[], char bmat, int n,
                  char* which, int nev, float tol, float resid[],
                  int ncv, float V[], int ldv, int iparam[],
                  int ipntr[], float workd[], float workl[],
                  int lworkl, int& info)

{

  int      irvec;
  logical* iselect;
  float*   iZ;

  irvec   = (int) rvec;
  iselect = new logical[ncv];
  iZ = (Z == NULL) ? &V[0] : Z;

  F77NAME(sneupd)(&irvec, &HowMny, iselect, dr, di, iZ, &ldz, &sigmar,
                  &sigmai, &workv[1], &bmat, &n, which, &nev, &tol,
                  resid, &ncv, &V[0], &ldv, &iparam[1], &ipntr[1],
                  &workd[1], &workl[1], &lworkl, &info HIDDEN_112);

  delete[] iselect;

} // neupp (float).  

inline void sauppError(int info)
{
  switch (info) {
  case     0:
    return;
  case    -8:
    throw ErrorExec("lapack error [sd]aupp ",info);
  case    -9:
    throw ErrorExec("arpack error START_RESID_ZERO  [sd]aupp ",info);
  case -9999:
    throw ErrorExec("arpack error ARNOLDI_NOT_BUILD  [sd]aupp ",info);
  case     1:
    throw ErrorExec("arpack error MAX_ITERATIONS  [sd]aupp ",info);
    return;
  case     3:
    throw ErrorExec("arpack error NO_SHIFTS_APPLIED  [sd]aupp ",info);
    return;
  default   :
    throw ErrorExec("arpack error AUPP_ERROR  [sd]aupp ",info);
  }

}

inline void seuppError(int info)
{

  switch (info) {
  case   0:
    return;
  case  -8:
  case  -9:
    throw ErrorExec("lapack LAPACK_ERROR [sd]eupp ",info);
  case -14:
    throw ErrorExec("lapack NOT_ACCURATE_EIG [sd]eupp ",info);
  case   1:
    throw ErrorExec("lapack REORDERING_ERROR [sd]eupp ",info);
  default :
    throw ErrorExec("lapack error [sd]eupp ",info);
  }

} // EuppError.


#endif // ARPACKFF_HPP
