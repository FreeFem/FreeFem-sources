

/*----------------------------------------------------------*/
/*															*/
/*						LPLIB	V3.00						*/
/*															*/
/*----------------------------------------------------------*/
/*															*/
/*	Description:		Handles threads, scheduling			*/
/*						& dependencies						*/
/*	Author:				Loic MARECHAL						*/
/*	Creation date:		feb 25 2008							*/
/*	Last modification:	jun 16 2010							*/
/*															*/
/*----------------------------------------------------------*/


/*----------------------------------------------------------*/
/* User available procedures' prototypes					*/
/*----------------------------------------------------------*/

int InitParallel(int);
void StopParallel(int);
int NewType(int, int);
void FreeType(int, int);
int BeginDependency(int, int, int);
void AddDependency(int, int, int);
void EndDependency(int, float [2]);
float LaunchParallel(int, int, int, void *, void *);
int LaunchPipeline(int, void *, void *, int, int *);
void WaitPipeline(int);
int ParallelMemClear(int , void *, long);
void ParallelQsort(int, void *, size_t, size_t, int (*)(const void *, const void *));
void HilbertRenumbering(int, int, double [6], double (*)[3], unsigned long long (*)[2]);
void HilbertRenumbering2D(int, int, double [4], double (*)[2], unsigned long long (*)[2]);
int AllocHash(int, int, int);
void FreeHash(int, int);
long long int AddHash(int, int, int, int, int, int, long long int);
