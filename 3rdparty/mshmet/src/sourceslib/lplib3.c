

/*----------------------------------------------------------*/
/*															*/
/*						LPLIB	V3.02						*/
/*															*/
/*----------------------------------------------------------*/
/*															*/
/*	Description:		Handles threads, scheduling			*/
/*						& dependencies						*/
/*	Author:				Loic MARECHAL						*/
/*	Creation date:		feb 25 2008							*/
/*	Last modification:	feb 15 2011							*/
/*															*/
/*----------------------------------------------------------*/


#ifndef SERIAL

/*----------------------------------------------------------*/
/* Includes													*/
/*----------------------------------------------------------*/

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <assert.h>
#include <errno.h>
#include "lplib3.h"
#ifdef __FreeBSD__
#include <sys/types.h>
#include <pmc.h>
#endif

/*----------------------------------------------------------*/
/* Defines													*/
/*----------------------------------------------------------*/

#define MaxLibPar 10
#define MaxPth 128
#define MaxTyp 100
#define DefNmbSmlBlk 64
#define DefNmbDepBlk 256
#define MaxTotPip 65536
#define MaxPipDep 100
#define MaxHsh 10

enum ParCmd {RunBigWrk, RunSmlWrk, ClrMem, EndPth};


/*----------------------------------------------------------*/
/* Structures' prototypes									*/
/*----------------------------------------------------------*/

typedef struct WrkSct
{
	int BegIdx, EndIdx, NmbDep, *DepWrdTab;
	int *DepIdxTab;
	struct WrkSct *pre, *nex;
}WrkSct;

typedef struct
{
	int NmbLin, NmbSmlWrk, NmbBigWrk, SmlWrkSiz, BigWrkSiz, DepWrkSiz, NmbDepWrd, *DepWrdMat, *DepIdxMat;
	char *RunDepTab;
	WrkSct *SmlWrkTab, *BigWrkTab;
}TypSct;

typedef struct
{
	int idx;
	char *ClrAdr;
	WrkSct *wrk;
	pthread_mutex_t mtx;
	pthread_cond_t cnd;
	pthread_t pth;
	struct ParSct *par;
}PthSct;

typedef struct PipSct
{
	int idx, NmbDep, DepTab[ MaxPipDep ];
	void *prc, *arg;
	pthread_t pth;
	struct ParSct *par;
}PipSct;

typedef struct BucSct
{
	int idx[3];
	long long int dat;
	struct BucSct *nex;
}BucSct;

typedef struct
{
	int mul[3], NmbOvf[ MaxPth ];
	TypSct *typ1, *typ2;
	BucSct *buc, *ovf[ MaxPth ];
}HshSct;

typedef struct ParSct
{
	int NmbCpu, WrkCpt, NmbPip, PenPip, RunPip, NmbTyp, BufMax, BufCpt, req, cmd, ClrLinSiz, *PipWrd;
	double sta[2];
	void (*prc)(int, int, int, void *), *arg;
	pthread_cond_t ParCnd, PipCnd;
	pthread_mutex_t ParMtx, PipMtx;
	pthread_t PipPth;
	PthSct *PthTab;
	TypSct *TypTab, *CurTyp, *DepTyp, *typ1, *typ2;
	WrkSct *NexWrk, *BufWrk[ MaxPth / 4 ];
	HshSct HshTab[ MaxHsh ];
}ParSct;

typedef struct
{
	unsigned long long (*idx)[2];
	double box[6], (*crd)[3], (*crd2)[2];
}ArgSct;


/*----------------------------------------------------------*/
/* Private procedures' prototypes							*/
/*----------------------------------------------------------*/

static int SetBit(int *, int);
static int GetBit(int *, int);
static int AndWrd(WrkSct *, char *);
static void SetWrd(WrkSct *, char *);
static void ClrWrd(WrkSct *, char *);
int CmpWrk(const void *, const void *);
static void *PipHdl(void *);
static void *PthHdl(void *);
static WrkSct *NexWrk(ParSct *, int);


/*----------------------------------------------------------*/
/* Global variables											*/
/*----------------------------------------------------------*/

ParSct *ParTab[ MaxLibPar+1 ];
int IniLibPar = 0;


/*----------------------------------------------------------*/
/* Init structures, scheduler and launch threads			*/
/*----------------------------------------------------------*/

int InitParallel(int NmbCpu)
{
	int i, ParIdx;
	ParSct *par = NULL;
	PthSct *pth;

	/* Check the number of requested cpu and clear the main par table at first call */

	if(NmbCpu > MaxPth)
		return(0);

	if(!IniLibPar)
	{
		IniLibPar = 1;

		for(i=1;i<=MaxLibPar;i++)
			ParTab[i] = NULL;
	}

	/* Allocate and build main parallel structure */

	for(ParIdx=1; ParIdx<=MaxLibPar; ParIdx++)
		if(!ParTab[ ParIdx ])
		{
			par = ParTab[ ParIdx ] = calloc(1, sizeof(ParSct));
			break;
		}

	if(!par)
		return(0);

	if(!(par->PthTab = calloc(NmbCpu, sizeof(PthSct))))
		return(0);

	if(!(par->TypTab = calloc((MaxTyp + 1), sizeof(TypSct))))
		return(0);

	if(!(par->PipWrd = calloc(MaxTotPip/32, sizeof(int))))
		return(0);

	par->NmbCpu = NmbCpu;
	par->WrkCpt = par->NmbPip = par->PenPip = par->RunPip = 0;

	/* Set the size of WP buffer */

	if(NmbCpu >= 4)
		par->BufMax = NmbCpu / 4;
	else
		par->BufMax = 1;

	pthread_mutex_init(&par->ParMtx, NULL);
	pthread_mutex_init(&par->PipMtx, NULL);
	pthread_cond_init(&par->ParCnd, NULL);
	pthread_cond_init(&par->PipCnd, NULL);

	/* Launch pthreads */

	for(i=0;i<par->NmbCpu;i++)
	{
		pth = &par->PthTab[i];
		pth->idx = i;
		pth->par = par;
		pthread_mutex_init(&pth->mtx, NULL);
		pthread_cond_init(&pth->cnd, NULL);
		pthread_create(&pth->pth, NULL, PthHdl, (void *)pth);
	}

	/* Wait for all threads to be up and wainting */

	pthread_mutex_lock(&par->ParMtx);

	while(par->WrkCpt < par->NmbCpu)
		pthread_cond_wait(&par->ParCnd, &par->ParMtx);

	pthread_mutex_unlock(&par->ParMtx);

	return(ParIdx);
}


/*----------------------------------------------------------*/
/* Stop all threads and free memories						*/
/*----------------------------------------------------------*/

void StopParallel(int ParIdx)
{
	int i;
	PthSct *pth;
	ParSct *par;

	/* Get and check lib parallel instance */

	if( (ParIdx < 1) || (ParIdx > MaxLibPar) || !(par = ParTab[ ParIdx ]) )
		return;

	/* Send stop to all threads */

	pthread_mutex_lock(&par->ParMtx);
	par->cmd = EndPth;
	pthread_mutex_unlock(&par->ParMtx);

	/* Wait for all threads to complete */

	for(i=0;i<par->NmbCpu;i++)
	{
		pth = &par->PthTab[i];
		pthread_mutex_lock(&pth->mtx);
		pthread_cond_signal(&pth->cnd);
		pthread_mutex_unlock(&pth->mtx);
		pthread_join(pth->pth, NULL);
	}

	pthread_mutex_destroy(&par->ParMtx);
	pthread_cond_destroy(&par->ParCnd);

	WaitPipeline(ParIdx);

	pthread_mutex_destroy(&par->PipMtx);
	pthread_cond_destroy(&par->PipCnd);

	/* Free memories */

	for(i=1;i<=MaxTyp;i++)
		if(par->TypTab[i].NmbLin)
			FreeType(ParIdx, i);

	free(par->PthTab);
	free(par->TypTab);
	free(par->PipWrd);
	free(par);

	ParTab[ ParIdx ] = NULL;
}


/*----------------------------------------------------------*/
/* Launch the loop prc on typ1 element depending on typ2	*/
/*----------------------------------------------------------*/

float LaunchParallel(int ParIdx, int TypIdx1, int TypIdx2, void *prc, void *PtrArg)
{
	int i;
	PthSct *pth;
	ParSct *par;
	TypSct *typ1, *typ2 = NULL;

	/* Get and check lib parallel instance */

	if( (ParIdx < 1) || (ParIdx > MaxLibPar) || !(par = ParTab[ ParIdx ]) )
		return(-1.);

	/* Check bounds */

	if( (TypIdx1 < 1) || (TypIdx1 > MaxTyp) || (TypIdx2 < 0) || (TypIdx2 > MaxTyp) || (TypIdx1 == TypIdx2) )
		return(-1.);

	typ1 =  &par->TypTab[ TypIdx1 ];

	if(TypIdx2)
	{
		/* Lock acces to global parameters */
    
		pthread_mutex_lock(&par->ParMtx);

		par->cmd = RunSmlWrk;
		par->prc = (void (*)(int, int, int, void *))prc;
		par->arg = PtrArg;
		par->typ1 = typ1;
		par->typ2 = typ2 = &par->TypTab[ TypIdx2 ];
		par->NexWrk = typ1->SmlWrkTab;
		par->BufCpt = 0;
		par->WrkCpt = 0;
		par->sta[0] = par->sta[1] = 0.;
    	par->req = 0;

		/* Clear running wp */

		for(i=0;i<par->NmbCpu;i++)
			par->PthTab[i].wrk = NULL;

		memset(typ1->RunDepTab, 0, typ1->NmbDepWrd * 32 * sizeof(char));

		/* Build a linked list of wp */

		for(i=0;i<par->typ1->NmbSmlWrk;i++)
		{
			typ1->SmlWrkTab[i].pre = &typ1->SmlWrkTab[ i-1 ];
			typ1->SmlWrkTab[i].nex = &typ1->SmlWrkTab[ i+1 ];
		}

		typ1->SmlWrkTab[0].pre = typ1->SmlWrkTab[ typ1->NmbSmlWrk - 1 ].nex = NULL;

		/* Start main loop : wake up threads and wait for completion or blocked threads */

		do
		{
			/* Search for some idle threads */

			par->req = 0;

			for(i=0;i<par->NmbCpu;i++)
			{
				pth = &par->PthTab[i];

				if(pth->wrk)
					continue;

				if(!(pth->wrk = NexWrk(par, i)))
				{
					par->req = 1;
					break;
				}

				/* Wake up the thread and provide it with a WP list */

				pthread_mutex_lock(&pth->mtx);
				pthread_cond_signal(&pth->cnd);
				pthread_mutex_unlock(&pth->mtx);
			}

			/* If every WP are done : exit the parallel loop */

			if(par->WrkCpt == typ1->NmbSmlWrk)
				break;

			/* Otherwise, wait for a blocked thread */

			pthread_cond_wait(&par->ParCnd, &par->ParMtx);
		}while(1);

		pthread_mutex_unlock(&par->ParMtx);

		/* Return the average speedup */

		return(par->sta[1] / par->sta[0]);
	}
	else
	{
		/* Lock acces to global parameters */

		pthread_mutex_lock(&par->ParMtx);

		par->cmd = RunBigWrk;
		par->prc = (void (*)(int, int, int, void *))prc;
		par->arg = PtrArg;
		par->typ1 = typ1;
		par->typ2 = NULL;
		par->WrkCpt = 0;

		for(i=0;i<typ1->NmbBigWrk;i++)
		{
			pth = &par->PthTab[i];
			pth->wrk = &typ1->BigWrkTab[i];
		}

		for(i=0;i<typ1->NmbBigWrk;i++)
		{
			pth = &par->PthTab[i];
			pthread_mutex_lock(&pth->mtx);
			pthread_cond_signal(&pth->cnd);
			pthread_mutex_unlock(&pth->mtx);
		}

		pthread_cond_wait(&par->ParCnd, &par->ParMtx);

		pthread_mutex_unlock(&par->ParMtx);

		/* Return the average speedup */

		return(par->NmbCpu);
	}
}


/*----------------------------------------------------------*/
/* Pthread handler, waits for job, does it, then signal end	*/
/*----------------------------------------------------------*/

static void *PthHdl(void *ptr)
{
	PthSct *pth = (PthSct *)ptr;
	ParSct *par = pth->par;

	/* Tell the scheduler if all threads are ready */

	pthread_mutex_lock(&par->ParMtx);
	par->WrkCpt++;
	pthread_cond_signal(&par->ParCnd);
	pthread_mutex_lock(&pth->mtx);
	pthread_mutex_unlock(&par->ParMtx);

	/* Enter main loop until StopParallel is send */

	do
	{
		/* Wait for a wake-up signal from the main loop */

		pthread_cond_wait(&pth->cnd, &pth->mtx);

		switch(par->cmd)
		{
			case RunBigWrk :
			{
				/* Launch a single big wp and signal completion to the scheduler */

				par->prc(pth->wrk->BegIdx, pth->wrk->EndIdx, pth->idx, par->arg);

				pthread_mutex_lock(&par->ParMtx);
				par->WrkCpt++;

				if(par->WrkCpt >= par->typ1->NmbBigWrk)
					pthread_cond_signal(&par->ParCnd);

				pthread_mutex_unlock(&par->ParMtx);
			}break;

			case RunSmlWrk :
			{
				do
				{
					/* Run the WP */

					par->prc(pth->wrk->BegIdx, pth->wrk->EndIdx, pth->idx, par->arg);

					/* Locked acces to global parameters : update WP count, tag WP done and signal the main loop */

					pthread_mutex_lock(&par->ParMtx);

					par->WrkCpt++;

					if(!(pth->wrk = NexWrk(par, pth->idx)))
					{
						par->req = 1;
						pthread_cond_signal(&par->ParCnd);
						pthread_mutex_unlock(&par->ParMtx);
						break;
					}

					if(par->req)
						pthread_cond_signal(&par->ParCnd);

					pthread_mutex_unlock(&par->ParMtx);
				}while(1);
			}break;

			case ClrMem :
			{
				/* Clear memory and signal completion to the scheduler */

				memset(pth->ClrAdr, 0, par->ClrLinSiz);

				pthread_mutex_lock(&par->ParMtx);
				par->WrkCpt++;
				pthread_cond_signal(&par->ParCnd);
				pthread_mutex_unlock(&par->ParMtx);
			}break;

			case EndPth :
			{
				/* Destroy the thread mutex and condition and call for join */

				pthread_mutex_unlock(&pth->mtx);
				pthread_mutex_destroy(&pth->mtx);
				pthread_cond_destroy(&pth->cnd);
				return(NULL);
			}break;
		}
	}while(1);

	return(NULL);
}


/*----------------------------------------------------------*/
/* Get the next WP to be computed							*/
/*----------------------------------------------------------*/

static WrkSct *NexWrk(ParSct *par, int PthIdx)
{
	int i;
	PthSct *pth = &par->PthTab[ PthIdx ];
	WrkSct *wrk;

	/* Update stats */

	par->sta[0]++;

	for(i=0;i<par->NmbCpu;i++)
		if(par->PthTab[i].wrk)
			par->sta[1]++;

	/* Remove previous work's tags */

	if(pth->wrk)
		ClrWrd(pth->wrk, par->typ1->RunDepTab);

	/* If the wp's buffer is empty search for some new compatible wp to fill in */

	if(!par->BufCpt)
	{
		wrk = par->NexWrk;

		while(wrk)
		{
			/* Check for dependencies */
    
			if(!AndWrd(wrk, par->typ1->RunDepTab))
			{
				par->BufWrk[ par->BufCpt++ ] = wrk;

				/* Unlink wp */

				if(wrk->pre)
					wrk->pre->nex = wrk->nex;
				else
					par->NexWrk = wrk->nex;

				if(wrk->nex)
					wrk->nex->pre = wrk->pre;

				/* Add new work's tags */

				SetWrd(wrk, par->typ1->RunDepTab);

				if(par->BufCpt == par->BufMax)
					break;
			}

			wrk = wrk->nex;
		}
	}

	/* Return the next available wp in buffer and unlink it from the todo list */

	if(par->BufCpt)
		return(par->BufWrk[ --par->BufCpt ]);
	else
		return(NULL);
}


/*----------------------------------------------------------*/
/* Allocate a new kind of elements and set work-packages	*/
/*----------------------------------------------------------*/

int NewType(int ParIdx, int NmbLin)
{
	int i, TypIdx=0, idx;
	TypSct *typ;
	ParSct *par;

	/* Get and check lib parallel instance */

	if( (ParIdx < 1) || (ParIdx > MaxLibPar) || !(par = ParTab[ ParIdx ]) )
		return(0);

	if(NmbLin < par->NmbCpu)
		return(0);

	/* Search for a free type structure */

	for(i=1;i<=MaxTyp;i++)
		if(!par->TypTab[i].NmbLin)
		{
			TypIdx = i;
			break;
		}

	if(!TypIdx)
		return(0);

	typ = &par->TypTab[ TypIdx ];
	typ->NmbLin = NmbLin;

	/* Compute the size of small work-packages */

	if(NmbLin >= DefNmbSmlBlk * par->NmbCpu)
		typ->SmlWrkSiz = NmbLin / (DefNmbSmlBlk * par->NmbCpu);
	else
		typ->SmlWrkSiz = NmbLin / par->NmbCpu;

	typ->NmbSmlWrk = NmbLin / typ->SmlWrkSiz;

	if(NmbLin != typ->NmbSmlWrk * typ->SmlWrkSiz)
		typ->NmbSmlWrk++;

	if(!(typ->SmlWrkTab = calloc(typ->NmbSmlWrk , sizeof(WrkSct))))
		return(0);

	/* Set small work-packages */

	idx = 0;

	for(i=0;i<typ->NmbSmlWrk;i++)
	{
		typ->SmlWrkTab[i].BegIdx = idx + 1;
		typ->SmlWrkTab[i].EndIdx = idx + typ->SmlWrkSiz;
		idx += typ->SmlWrkSiz;
	}

	typ->SmlWrkTab[ typ->NmbSmlWrk - 1 ].EndIdx = NmbLin;

	/* Compute the size of big work-packages */

	typ->BigWrkSiz = NmbLin / par->NmbCpu;
	typ->NmbBigWrk = par->NmbCpu;

	if(!(typ->BigWrkTab = calloc(typ->NmbBigWrk , sizeof(WrkSct))))
		return(0);

	/* Set big work-packages */

	idx = 0;

	for(i=0;i<typ->NmbBigWrk;i++)
	{
		typ->BigWrkTab[i].BegIdx = idx + 1;
		typ->BigWrkTab[i].EndIdx = idx + typ->BigWrkSiz;
		idx += typ->BigWrkSiz;
	}

	typ->BigWrkTab[ typ->NmbBigWrk - 1 ].EndIdx = NmbLin;

	return(TypIdx);
}


/*----------------------------------------------------------*/
/* Add this kind of element to the free-list				*/
/*----------------------------------------------------------*/

void FreeType(int ParIdx, int TypIdx)
{
	TypSct *typ;
	ParSct *par;

	/* Get and check lib parallel instance */

	if( (ParIdx < 1) || (ParIdx > MaxLibPar) || !(par = ParTab[ ParIdx ]) )
		return;

	/* Check bounds and free mem */

	if( (TypIdx < 1) || (TypIdx > MaxTyp) )
		return;

	typ = &par->TypTab[ TypIdx ];

	if(typ->SmlWrkTab)
		free(typ->SmlWrkTab);

	if(typ->BigWrkTab)
		free(typ->BigWrkTab);

	if(typ->DepIdxMat)
		free(typ->DepIdxMat);

	if(typ->RunDepTab)
		free(typ->RunDepTab);

	if(typ->DepWrdMat)
		free(typ->DepWrdMat);

	memset(typ, 0, sizeof(TypSct));
}


/*----------------------------------------------------------*/
/* Allocate a dependency matrix linking both types			*/
/*----------------------------------------------------------*/

int BeginDependency(int ParIdx, int TypIdx1, int TypIdx2)
{
	int i;
	ParSct *par;
	TypSct *typ1, *typ2;

	/* Get and check lib parallel instance */

	if( (ParIdx < 1) || (ParIdx > MaxLibPar) || !(par = ParTab[ ParIdx ]) )
		return(0);

	/* Check bounds */

	par->CurTyp = typ1 = &par->TypTab[ TypIdx1 ];
	par->DepTyp = typ2 = &par->TypTab[ TypIdx2 ];

	if( (TypIdx1 < 1) || (TypIdx1 > MaxTyp) || (TypIdx2 < 1) || (TypIdx2 > MaxTyp) || (typ1 == typ2) \
	|| !typ1->NmbLin || !typ2->NmbLin)
	{
		return(0);
	}

	/* Compute dependency table's size */

	if(typ2->NmbLin >= DefNmbDepBlk * par->NmbCpu)
		typ1->DepWrkSiz = typ2->NmbLin / (DefNmbDepBlk * par->NmbCpu);
	else
		typ1->DepWrkSiz = typ2->NmbLin / par->NmbCpu;

	typ1->NmbDepWrd = typ2->NmbLin / (typ1->DepWrkSiz * 32);

	if(typ2->NmbLin != typ1->NmbDepWrd * typ1->DepWrkSiz * 32)
		typ1->NmbDepWrd++;

	/* Allocate a global dependency table */

	if(!(typ1->DepWrdMat = calloc(typ1->NmbSmlWrk * typ1->NmbDepWrd, sizeof(int))))
		return(0);

	/* Then spread sub-tables among WP */

	for(i=0;i<typ1->NmbSmlWrk;i++)
	{
		typ1->SmlWrkTab[i].NmbDep = 0;
		typ1->SmlWrkTab[i].DepWrdTab = &typ1->DepWrdMat[ i * typ1->NmbDepWrd ];
	}

	/* Allocate a running tags table */

	if(!(typ1->RunDepTab = calloc(typ1->NmbDepWrd * 32, sizeof(char))))
		return(0);

	return(1);
}


/*----------------------------------------------------------*/
/* Type1 element idx1 depends on type2 element idx2			*/
/*----------------------------------------------------------*/

void AddDependency(int ParIdx, int idx1, int idx2)
{
	WrkSct *wrk;
	ParSct *par;

	/* Get and check lib parallel instance */

	if( (ParIdx < 1) || (ParIdx > MaxLibPar) || !(par = ParTab[ ParIdx ]) )
		return;

	/* Set and count dependency bit */

	wrk = &par->CurTyp->SmlWrkTab[ (idx1-1) / par->CurTyp->SmlWrkSiz ];

	if(!SetBit(wrk->DepWrdTab, (idx2-1) / par->CurTyp->DepWrkSiz ))
		wrk->NmbDep++;
}


/*----------------------------------------------------------*/
/* Sort wp depending on their number of dependencies		*/
/*----------------------------------------------------------*/

void EndDependency(int ParIdx, float DepSta[2])
{
	int i, j, idx=0, NmbDep, NmbDepBit, TotNmbDep = 0;
	ParSct *par;
	TypSct *typ1, *typ2;
	WrkSct *wrk;

	/* Get and check lib parallel instance */

	if( (ParIdx < 1) || (ParIdx > MaxLibPar) || !(par = ParTab[ ParIdx ]) )
		return;

	/* Compute average number of collisions */

	DepSta[1] = 0.;
	typ1 = par->CurTyp;
	typ2 = par->DepTyp;

	for(i=0;i<typ1->NmbSmlWrk;i++)
	{
		TotNmbDep += typ1->SmlWrkTab[i].NmbDep;

		if(typ1->SmlWrkTab[i].NmbDep > DepSta[1])
			DepSta[1] = typ1->SmlWrkTab[i].NmbDep;
	}

	DepSta[0] = TotNmbDep;

	/* Allocate a global dependency index table */

	if(!(typ1->DepIdxMat = calloc(TotNmbDep, sizeof(int))))
		return;

	/* Then spread and fill the sub-tables among WP */

	NmbDep = typ1->NmbDepWrd * 32;

	for(i=0;i<typ1->NmbSmlWrk;i++)
	{
		wrk = &typ1->SmlWrkTab[i];
		wrk->DepIdxTab = &typ1->DepIdxMat[ idx ];
		idx += wrk->NmbDep;
		wrk->NmbDep = 0;

		for(j=0;j<NmbDep;j++)
			if(GetBit(wrk->DepWrdTab, j))
				wrk->DepIdxTab[ wrk->NmbDep++ ] = j;
	}

	/* Compute stats */

	NmbDepBit = typ2->NmbLin / typ1->DepWrkSiz;

	if(typ2->NmbLin - NmbDepBit * typ1->DepWrkSiz)
		NmbDepBit++;

	DepSta[0] = 100 * DepSta[0] / (typ1->NmbSmlWrk * NmbDepBit);
	DepSta[1] = 100 * DepSta[1] / NmbDepBit;

	/* Sort WP from highest collision number to the lowest */

	qsort(typ1->SmlWrkTab, typ1->NmbSmlWrk, sizeof(WrkSct), CmpWrk);

	par->CurTyp = NULL;
}


/*----------------------------------------------------------*/
/* Test and set a bit in a multibyte number					*/
/*----------------------------------------------------------*/

static int SetBit(int *tab, int idx)
{
	int res = ( tab[ idx >> 5 ] & (1 << (idx & 31)) );
	tab[ idx >> 5 ] |= 1 << (idx & 31);
	return(res);
}


/*----------------------------------------------------------*/
/* Test a bit in a multibyte number							*/
/*----------------------------------------------------------*/

static int GetBit(int *tab, int idx)
{
	return( tab[ idx >> 5 ] & (1 << (idx & 31)) );
}


/*----------------------------------------------------------*/
/* Check wether two WP share common resources -> locked		*/
/*----------------------------------------------------------*/

static int AndWrd(WrkSct *wrk, char *wrd)
{
	int i;

	for(i=0;i<wrk->NmbDep;i++)
		if(wrd[ wrk->DepIdxTab[i] ])
			return(1);

	return(0);
}

static void SetWrd(WrkSct *wrk, char *wrd)
{
	int i;

	for(i=0;i<wrk->NmbDep;i++)
		wrd[ wrk->DepIdxTab[i] ] = 1;
}

static void ClrWrd(WrkSct *wrk, char *wrd)
{
	int i;

	for(i=0;i<wrk->NmbDep;i++)
		wrd[ wrk->DepIdxTab[i] ] = 0;
}


/*----------------------------------------------------------*/
/* Compare two workpackages number of bits					*/
/*----------------------------------------------------------*/

int CmpWrk(const void *ptr1, const void *ptr2)
{
	WrkSct *w1, *w2;

	w1 = (WrkSct *)ptr1;
	w2 = (WrkSct *)ptr2;

	if(w1->NmbDep > w2->NmbDep)
		return(-1);
	else if(w1->NmbDep < w2->NmbDep)
		return(1);
	else
		return(0);
}


/*----------------------------------------------------------*/
/* Launch the loop prc on typ1 element depending on typ2	*/
/*----------------------------------------------------------*/

int ParallelMemClear(int ParIdx, void *PtrArg, long siz)
{
	char *tab = (char *)PtrArg;
	int i;
	PthSct *pth;
	ParSct *par;

	/* Get and check lib parallel instance, adresse and size */

	if( (ParIdx < 1) || (ParIdx > MaxLibPar) || !(par = ParTab[ ParIdx ]) || !tab || (siz < par->NmbCpu) )
		return(0);

	/* Lock acces to global parameters */
	
	pthread_mutex_lock(&par->ParMtx);

	par->cmd = ClrMem;
	par->ClrLinSiz = siz / par->NmbCpu;
	par->WrkCpt = 0;

	/* Spread the buffer among each thread and wake then up */

	for(i=0;i<par->NmbCpu;i++)
	{
		pth = &par->PthTab[i];
		pth->ClrAdr = &tab[ i * par->ClrLinSiz ];

		pthread_mutex_lock(&pth->mtx);
		pthread_cond_signal(&pth->cnd);
		pthread_mutex_unlock(&pth->mtx);
	}

	/* Wait for each thread to complete */

	while(par->WrkCpt < par->NmbCpu)
		pthread_cond_wait(&par->ParCnd, &par->ParMtx);

	pthread_mutex_unlock(&par->ParMtx);

	return(1);
}


/*----------------------------------------------------------*/
/* Wait for a condition, launch and detach a user procedure	*/
/*----------------------------------------------------------*/

int LaunchPipeline(int ParIdx, void *prc, void *PtrArg, int NmbDep, int *DepTab)
{
	int i;
	PipSct *NewPip=NULL;
	ParSct *par;

	/* Get and check lib parallel instance and the number of pipes and dependencies */

	if( (ParIdx < 1) || (ParIdx > MaxLibPar) || !(par = ParTab[ ParIdx ]) \
	||	(NmbDep > MaxPipDep) || (par->NmbPip >= MaxTotPip) )
	{
		return(0);
	}

	/* Allocate and setup a new pipe */

	if(!(NewPip = malloc(sizeof(PipSct))))
		return(0);

	NewPip->prc = prc;
	NewPip->arg = PtrArg;
	NewPip->par = par;
	NewPip->NmbDep = NmbDep;

	for(i=0;i<NmbDep;i++)
		NewPip->DepTab[i] = DepTab[i];

	/* Lock pipe mutex, increment pipe counter and launch the pipe regardless dependencies */

	pthread_mutex_lock(&par->PipMtx);
	NewPip->idx = ++par->NmbPip;
	par->PenPip++;
	pthread_create(&NewPip->pth, NULL, PipHdl, (void *)NewPip);
	pthread_mutex_unlock(&par->PipMtx);

	return(NewPip->idx);
}


/*----------------------------------------------------------*/
/* Thread handler launching and waitint for user's procedure*/
/*----------------------------------------------------------*/

static void *PipHdl(void *ptr)
{
	int RunFlg, i;
	PipSct *pip = (PipSct *)ptr;
	ParSct *par = pip->par;
	void (*prc)(void *);

	/* Wait for conditions to be met */

	do
	{
		pthread_mutex_lock(&par->PipMtx);

		if(par->RunPip < par->NmbCpu)
		{
			RunFlg = 1;

			for(i=0;i<pip->NmbDep;i++)
				if(!GetBit(par->PipWrd, pip->DepTab[i]))
				{
					RunFlg = 0;
					break;
				}
		}

		if(!RunFlg)
		{
			pthread_mutex_unlock(&par->PipMtx);
			usleep(1000);
		}
	}while(!RunFlg);

	/* Execute the user's procedure and set the flag to 2 (done) */

	prc = (void (*)(void *))pip->prc;
	par->RunPip++;

	pthread_mutex_unlock(&par->PipMtx);

	prc(pip->arg);

	pthread_mutex_lock(&par->PipMtx);
	SetBit(par->PipWrd, pip->idx);
	par->PenPip--;
	par->RunPip--;
	free(pip);
	pthread_mutex_unlock(&par->PipMtx);

	return(NULL);
}


/*----------------------------------------------------------*/
/* Wait for all pipelined procedures to complete			*/
/*----------------------------------------------------------*/

void WaitPipeline(int ParIdx)
{
	int PenPip;
	ParSct *par;

	/* Get and check lib parallel instance */

	if( (ParIdx < 1) || (ParIdx > MaxLibPar) || !(par = ParTab[ ParIdx ]) )
		return;

	do
	{
		pthread_mutex_lock(&par->PipMtx);
		PenPip = par->PenPip;
		pthread_mutex_unlock(&par->PipMtx);
		usleep(1000);
	}while(PenPip);
}


/*-
 * Copyright (c) 1992, 1993
 *	The Regents of the University of California.  All rights reserved.
 * Multithread implementation Copyright (c) 2006, 2007 Diomidis Spinellis.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 4. Neither the name of the University nor the names of its contributors
 *    may be used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 */

#define verify(x) (x)
#define DLOG(...)
#define min(a, b)	(a) < (b) ? a : b

typedef int		 cmp_t(const void *, const void *);
static inline char	*med3(char *, char *, char *, cmp_t *, void *);
static inline void	 swapfunc(char *, char *, int, int);

/*
 * Qsort routine from Bentley & McIlroy's "Engineering a Sort Function".
 */
#define swapcode(TYPE, parmi, parmj, n) { 		\
	long i = (n) / sizeof (TYPE); 			\
	TYPE *pi = (TYPE *) (parmi); 		\
	TYPE *pj = (TYPE *) (parmj); 		\
	do { 						\
		TYPE	t = *pi;		\
		*pi++ = *pj;				\
		*pj++ = t;				\
        } while (--i > 0);				\
}


static inline void swapfunc(char *a, char *b, int n, int swaptype)
{
	if(swaptype <= 1)
		swapcode(long, a, b, n)
	else
		swapcode(char, a, b, n)
}

#define swap(a, b)					\
	if (swaptype == 0) {				\
		long t = *(long *)(a);			\
		*(long *)(a) = *(long *)(b);		\
		*(long *)(b) = t;			\
	} else						\
		swapfunc(a, b, es, swaptype)

#define vecswap(a, b, n) 	if ((n) > 0) swapfunc(a, b, n, swaptype)

#define	CMP(t, x, y) (cmp((x), (y)))

static inline char *med3(char *a, char *b, char *c, cmp_t *cmp, void *thunk)
{
	return CMP(thunk, a, b) < 0 ?
	       (CMP(thunk, b, c) < 0 ? b : (CMP(thunk, a, c) < 0 ? c : a ))
              :(CMP(thunk, b, c) > 0 ? b : (CMP(thunk, a, c) < 0 ? a : c ));
}

/*
 * We use some elaborate condition variables and signalling
 * to ensure a bound of the number of active threads at
 * 2 * maxthreads and the size of the thread data structure
 * to maxthreads.
 */

/* Condition of starting a new thread. */
enum thread_state {
	ts_idle,		/* Idle, waiting for instructions. */
	ts_work,		/* Has work to do. */
	ts_term			/* Asked to terminate. */
};

/* Variant part passed to qsort invocations. */
struct qsort {
	enum thread_state st;	/* For coordinating work. */
	struct common *common;	/* Common shared elements. */
	void *a;		/* Array base. */
	size_t n;		/* Number of elements. */
	pthread_t id;		/* Thread id. */
	pthread_mutex_t mtx_st;	/* For signalling state change. */
	pthread_cond_t cond_st;	/* For signalling state change. */
};

/* Invariant common part, shared across invocations. */
struct common {
	int swaptype;		/* Code to use for swapping */
	size_t es;		/* Element size. */
	void *thunk;		/* Thunk for qsort_r */
	cmp_t *cmp;		/* Comparison function */
	int nthreads;		/* Total number of pool threads. */
	int idlethreads;	/* Number of idle threads in pool. */
	int forkelem;		/* Minimum number of elements for a new thread. */
	struct qsort *pool;	/* Fixed pool of threads. */
	pthread_mutex_t mtx_al;	/* For allocating threads in the pool. */
};

static void *qsort_thread(void *p);

/* The multithreaded qsort public interface */
void qsort_mt(void *a, size_t n, size_t es, cmp_t *cmp, int maxthreads, int forkelem)
{
	struct qsort *qs;
	struct common c;
	int i, islot;
	bool bailout = true;

	if (n < forkelem)
		goto f1;
	errno = 0;
#ifdef __FreeBSD__
	if (maxthreads == 0) {
		/*
		 * Other candidates:
		 * NPROC environment variable (BSD/OS, CrayOS)
		 * sysctl hw.ncpu or kern.smp.cpus
		 */
	        int ncpu; 
		if (pmc_init() == 0 && (ncpu = pmc_ncpu()) != -1)
			maxthreads = ncpu;
		else
			maxthreads = 2;
	}
#endif
	/* XXX temporarily disabled for stress and performance testing.
	if (maxthreads == 1)
		goto f1;
	*/
	/* Try to initialize the resources we need. */
	if (pthread_mutex_init(&c.mtx_al, NULL) != 0)
		goto f1;
	if ((c.pool = (struct qsort *)calloc(maxthreads, sizeof(struct qsort))) ==NULL)
		goto f2;
	for (islot = 0; islot < maxthreads; islot++) {
		qs = &c.pool[islot];
		if (pthread_mutex_init(&qs->mtx_st, NULL) != 0)
			goto f3;
		if (pthread_cond_init(&qs->cond_st, NULL) != 0) {
			verify(pthread_mutex_destroy(&qs->mtx_st));
			goto f3;
		}
		qs->st = ts_idle;
		qs->common = &c;
		if (pthread_create(&qs->id, NULL, qsort_thread, qs) != 0) {
			verify(pthread_mutex_destroy(&qs->mtx_st));
			verify(pthread_cond_destroy(&qs->cond_st));
			goto f3;
		}
	}

	/* All systems go. */
	bailout = false;

	/* Initialize common elements. */
	c.swaptype = ((char *)a - (char *)0) % sizeof(long) || \
		es % sizeof(long) ? 2 : es == sizeof(long)? 0 : 1;
	c.es = es;
	c.cmp = cmp;
	c.forkelem = forkelem;
	c.idlethreads = c.nthreads = maxthreads;

	/* Hand out the first work batch. */
	qs = &c.pool[0];
	verify(pthread_mutex_lock(&qs->mtx_st));
	qs->a = a;
	qs->n = n;
	qs->st = ts_work;
	c.idlethreads--;
	verify(pthread_cond_signal(&qs->cond_st));
	verify(pthread_mutex_unlock(&qs->mtx_st));

	/* 
	 * Wait for all threads to finish, and
	 * free acquired resources.
	 */
f3:	for (i = 0; i < islot; i++) {
		qs = &c.pool[i];
		if (bailout) {
			verify(pthread_mutex_lock(&qs->mtx_st));
			qs->st = ts_term;
			verify(pthread_cond_signal(&qs->cond_st));
			verify(pthread_mutex_unlock(&qs->mtx_st));
		}
		verify(pthread_join(qs->id, NULL));
		verify(pthread_mutex_destroy(&qs->mtx_st));
		verify(pthread_cond_destroy(&qs->cond_st));
	}
	free(c.pool);
f2:	verify(pthread_mutex_destroy(&c.mtx_al));
	if (bailout) {
		DLOG("Resource initialization failed; bailing out.\n");
		/* XXX should include a syslog call here */
		fprintf(stderr, "Resource initialization failed; bailing out.\n");
f1:		qsort(a, n, es, cmp);
	}
}

#define thunk NULL

/*
 * Allocate an idle thread from the pool, lock its
 * mutex, change its state to work, decrease the number
 * of idle threads, and return a
 * pointer to its data area.
 * Return NULL, if no thread is available.
 */
static struct qsort *allocate_thread(struct common *c)
{
	int i;

	verify(pthread_mutex_lock(&c->mtx_al));
	for (i = 0; i < c->nthreads; i++)
		if (c->pool[i].st == ts_idle) {
			c->idlethreads--;
			c->pool[i].st = ts_work;
			verify(pthread_mutex_lock(&c->pool[i].mtx_st));
			verify(pthread_mutex_unlock(&c->mtx_al));
			return (&c->pool[i]);
		}
	verify(pthread_mutex_unlock(&c->mtx_al));
	return (NULL);
}

/* Thread-callable quicksort. */
static void qsort_algo(struct qsort *qs)
{
	char *pa, *pb, *pc, *pd, *pl, *pm, *pn;
	int d, r, swaptype, swap_cnt;
	void *a;			/* Array of elements. */
	size_t n, es;			/* Number of elements; size. */
	cmp_t *cmp;
	int nl, nr;
	struct common *c;
	struct qsort *qs2;
	pthread_t id;

	/* Initialize qsort arguments. */
	id = qs->id;
	c = qs->common;
	es = c->es;
	cmp = c->cmp;
	swaptype = c->swaptype;
	a = qs->a;
	n = qs->n;
top:
	DLOG("%10x n=%-10d Sort starting.\n", id, n);

	/* From here on qsort(3) business as usual. */
	swap_cnt = 0;
	if (n < 7) {
		for (pm = (char *)a + es; pm < (char *)a + n * es; pm += es)
			for (pl = pm;
			     pl > (char *)a && CMP(thunk, pl - es, pl) > 0;
			     pl -= es)
				swap(pl, pl - es);
		return;
	}
	pm = (char *)a + (n / 2) * es;
	if (n > 7) {
		pl = a;
		pn = (char *)a + (n - 1) * es;
		if (n > 40) {
			d = (n / 8) * es;
			pl = med3(pl, pl + d, pl + 2 * d, cmp, thunk);
			pm = med3(pm - d, pm, pm + d, cmp, thunk);
			pn = med3(pn - 2 * d, pn - d, pn, cmp, thunk);
		}
		pm = med3(pl, pm, pn, cmp, thunk);
	}
	swap(a, pm);
	pa = pb = (char *)a + es;

	pc = pd = (char *)a + (n - 1) * es;
	for (;;) {
		while (pb <= pc && (r = CMP(thunk, pb, a)) <= 0) {
			if (r == 0) {
				swap_cnt = 1;
				swap(pa, pb);
				pa += es;
			}
			pb += es;
		}
		while (pb <= pc && (r = CMP(thunk, pc, a)) >= 0) {
			if (r == 0) {
				swap_cnt = 1;
				swap(pc, pd);
				pd -= es;
			}
			pc -= es;
		}
		if (pb > pc)
			break;
		swap(pb, pc);
		swap_cnt = 1;
		pb += es;
		pc -= es;
	}
	if (swap_cnt == 0) {  /* Switch to insertion sort */
		for (pm = (char *)a + es; pm < (char *)a + n * es; pm += es)
			for (pl = pm;
			     pl > (char *)a && CMP(thunk, pl - es, pl) > 0;
			     pl -= es)
				swap(pl, pl - es);
		return;
	}

	pn = (char *)a + n * es;
	r = min(pa - (char *)a, pb - pa);
	vecswap(a, pb - r, r);
	r = min(pd - pc, pn - pd - es);
	vecswap(pb, pn - r, r);

	nl = (pb - pa) / es;
	nr = (pd - pc) / es;
	DLOG("%10x n=%-10d Partitioning finished ln=%d rn=%d.\n", id, n, nl, nr);

	/* Now try to launch subthreads. */
	if (nl > c->forkelem && nr > c->forkelem &&
	    (qs2 = allocate_thread(c)) != NULL) {
		DLOG("%10x n=%-10d Left farmed out to %x.\n", id, n, qs2->id);
		qs2->a = a;
		qs2->n = nl;
		verify(pthread_cond_signal(&qs2->cond_st));
		verify(pthread_mutex_unlock(&qs2->mtx_st));
	} else if (nl > 0) {
		DLOG("%10x n=%-10d Left will be done in-house.\n", id, n);
		qs->a = a;
		qs->n = nl;
		qsort_algo(qs);
	}
	if (nr > 0) {
		DLOG("%10x n=%-10d Right will be done in-house.\n", id, n);
		a = pn - nr * es;
		n = nr;
		goto top;
	}
}

/* Thread-callable quicksort. */
static void *qsort_thread(void *p)
{
	struct qsort *qs, *qs2;
	int i;
	struct common *c;
	pthread_t id;

	qs = p;
	id = qs->id;
	c = qs->common;
again:
	/* Wait for work to be allocated. */
	DLOG("%10x n=%-10d Thread waiting for work.\n", id, 0);
	verify(pthread_mutex_lock(&qs->mtx_st));
	while (qs->st == ts_idle)
		verify(pthread_cond_wait(&qs->cond_st, &qs->mtx_st));
	verify(pthread_mutex_unlock(&qs->mtx_st));
	if (qs->st == ts_term) {
		DLOG("%10x n=%-10d Thread signalled to exit.\n", id, 0);
		return(NULL);
	}
	assert(qs->st == ts_work);

	qsort_algo(qs);

	verify(pthread_mutex_lock(&c->mtx_al));
	qs->st = ts_idle;
	c->idlethreads++;
	DLOG("%10x n=%-10d Finished idlethreads=%d.\n", id, 0, c->idlethreads);
	if (c->idlethreads == c->nthreads) {
		DLOG("%10x n=%-10d All threads idle, signalling shutdown.\n", id, 0);
		for (i = 0; i < c->nthreads; i++) {
			qs2 = &c->pool[i];
			if (qs2 == qs)
				continue;
			verify(pthread_mutex_lock(&qs2->mtx_st));
			qs2->st = ts_term;
			verify(pthread_cond_signal(&qs2->cond_st));
			verify(pthread_mutex_unlock(&qs2->mtx_st));
		}
		DLOG("%10x n=%-10d Shutdown signalling complete.\n", id, 0);
		verify(pthread_mutex_unlock(&c->mtx_al));
		return(NULL);
	}
	verify(pthread_mutex_unlock(&c->mtx_al));
	goto again;
}


/*----------------------------------------------------------*/
/* Lplib qsort_mt encapsulation								*/
/*----------------------------------------------------------*/

void ParallelQsort(int ParIdx, void *base, size_t nel, size_t width, int (*compar)(const void *, const void *))
{
	ParSct *par;

	/* Get and check lib parallel instance */

	if( (ParIdx < 1) || (ParIdx > MaxLibPar) || !(par = ParTab[ ParIdx ]) )
		return;

	qsort_mt(base, nel, width, compar, par->NmbCpu, 10000);
}


/*----------------------------------------------------------*/
/* Compute the hilbert code from 3d coordinates				*/
/*----------------------------------------------------------*/

static void RenPrc(int BegIdx, int EndIdx, int PthIdx, ArgSct *arg)
{
	unsigned long long IntCrd[3], m=1LL<<63, cod;
	int i, j, b, GeoWrd, NewWrd, BitTab[3] = {1,2,4};
	double dbl;
	int rot[8], GeoCod[8]={0,3,7,4,1,2,6,5};  /* Z curve = {5,4,7,6,1,0,3,2} */
	int HilCod[8][8] = {{0,7,6,1,2,5,4,3}, {0,3,4,7,6,5,2,1}, {0,3,4,7,6,5,2,1}, {2,3,0,1,6,7,4,5},\
						{2,3,0,1,6,7,4,5}, {6,5,2,1,0,3,4,7}, {6,5,2,1,0,3,4,7}, {4,3,2,5,6,1,0,7}};

	for(i=BegIdx; i<=EndIdx; i++)
	{
		/* Convert double precision coordinates to integers */

		for(j=0;j<3;j++)
		{
			dbl = (arg->crd[i][j] - arg->box[j]) * arg->box[j+3];
			IntCrd[j] = dbl;
		}

		/* Binary hilbert renumbering loop */

		cod = 0;

		for(j=0;j<8;j++)
			rot[j] = GeoCod[j];

		for(b=0;b<21;b++)
		{
			GeoWrd = 0;

			for(j=0;j<3;j++)
			{
				if(IntCrd[j] & m)
					GeoWrd |= BitTab[j];

				IntCrd[j] = IntCrd[j]<<1;
			}

			NewWrd = rot[ GeoWrd ];
			cod = cod<<3 | NewWrd;

			for(j=0;j<8;j++)
				rot[j] = HilCod[ NewWrd ][ rot[j] ];
		}

		arg->idx[i][0] = cod;
		arg->idx[i][1] = i;
	}
}


/*----------------------------------------------------------*/
/* Comparison of two items for the qsort					*/
/*----------------------------------------------------------*/

int CmpPrc(const void *a, const void *b)
{
	unsigned long long *pa = (unsigned long long *)a, *pb = (unsigned long long *)b;

	if(*pa > *pb)
		return(1);
	else
		return(-1);
}


/*----------------------------------------------------------*/
/* Renumber a set of coordinates through a Hilbert SFC		*/
/*----------------------------------------------------------*/

void HilbertRenumbering(int ParIdx, int NmbLin, double box[6], double (*crd)[3], unsigned long long (*idx)[2])
{
	int i, NewTyp;
	double len = pow(2,64);
	ParSct *par;
	ArgSct arg;

	/* Get and check lib parallel instance */

	if( (ParIdx < 1) || (ParIdx > MaxLibPar) || !(par = ParTab[ ParIdx ]) )
		return;

	NewTyp = NewType(ParIdx, NmbLin);
	arg.crd = crd;
	arg.idx = idx;
	arg.box[0] = box[0];
	arg.box[1] = box[1];
	arg.box[2] = box[2];
	arg.box[3] = len / (box[3] - box[0]);
	arg.box[4] = len / (box[4] - box[1]);
	arg.box[5] = len / (box[5] - box[2]);

	LaunchParallel(ParIdx, NewTyp, 0, (void *)RenPrc, (void *)&arg);
	ParallelQsort(ParIdx, &idx[1][0], NmbLin, 2 * sizeof(long long), CmpPrc);

	for(i=1;i<=NmbLin;i++)
		idx[ idx[i][1] ][0] = i;
}


/*----------------------------------------------------------*/
/* Compute the hilbert code from 2d coordinates				*/
/*----------------------------------------------------------*/

static void RenPrc2D(int BegIdx, int EndIdx, int PthIdx, ArgSct *arg)
{
	unsigned long long IntCrd[2], m=1LL<<62, cod;
	int i, j, b, GeoWrd, NewWrd, BitTab[2] = {1,2};
	double dbl;
	int rot[4], GeoCod[4]={1,2,0,3};
	int HilCod[4][4] = {{0,3,2,1}, {0,1,2,3}, {0,1,2,3}, {2,1,0,3}};

	for(i=BegIdx; i<=EndIdx; i++)
	{
		/* Convert double precision coordinates to integers */

		for(j=0;j<2;j++)
		{
			dbl = (arg->crd2[i][j] - arg->box[j]) * arg->box[j+2];
			IntCrd[j] = dbl;
		}

		/* Binary hilbert renumbering loop */

		cod = 0;

		for(j=0;j<4;j++)
			rot[j] = GeoCod[j];

		for(b=0;b<31;b++)
		{
			GeoWrd = 0;

			for(j=0;j<2;j++)
			{
				if(IntCrd[j] & m)
					GeoWrd |= BitTab[j];

				IntCrd[j] = IntCrd[j]<<1;
			}

			NewWrd = rot[ GeoWrd ];
			cod = cod<<2 | NewWrd;

			for(j=0;j<4;j++)
				rot[j] = HilCod[ NewWrd ][ rot[j] ];
		}

		arg->idx[i][0] = cod;
		arg->idx[i][1] = i;
	}
}


/*----------------------------------------------------------*/
/* Renumber a set of 2D coordinates through a Hilbert SFC	*/
/*----------------------------------------------------------*/

void HilbertRenumbering2D(int ParIdx, int NmbLin, double box[4], double (*crd)[2], unsigned long long (*idx)[2])
{
	int i, NewTyp;
	double len = pow(2,62);
	ParSct *par;
	ArgSct arg;

	/* Get and check lib parallel instance */

	if( (ParIdx < 1) || (ParIdx > MaxLibPar) || !(par = ParTab[ ParIdx ]) )
		return;

	NewTyp = NewType(ParIdx, NmbLin);
	arg.crd2 = crd;
	arg.idx = idx;
	arg.box[0] = box[0];
	arg.box[1] = box[1];
	arg.box[2] = len / (box[2] - box[0]);
	arg.box[3] = len / (box[3] - box[1]);

	LaunchParallel(ParIdx, NewTyp, 0, (void *)RenPrc2D, (void *)&arg);
	ParallelQsort(ParIdx, &idx[1][0], NmbLin, 2 * sizeof(long long), CmpPrc);

	for(i=1;i<=NmbLin;i++)
		idx[ idx[i][1] ][0] = i;
}


/*----------------------------------------------------------*/
/* Setup a parallel hash table and return its index			*/
/*----------------------------------------------------------*/

int AllocHash(int ParIdx, int BasTyp, int DepTyp)
{
	int HshIdx, mul, i;
	ParSct *par;
	HshSct *hsh;
	TypSct *typ1, *typ2;

	/* Get and check lib parallel instance */

	if( (ParIdx < 1) || (ParIdx > MaxLibPar) || !(par = ParTab[ ParIdx ]) || !BasTyp || !DepTyp)
		return(0);

	typ1 = &par->TypTab[ BasTyp ];
	typ2 = &par->TypTab[ DepTyp ];

	for(HshIdx=1; HshIdx< MaxHsh; HshIdx++)
		if(!par->HshTab[ HshIdx ].typ1)
			break;

	if( (HshIdx == MaxHsh) || !typ1->NmbLin || !typ2->NmbLin )
		return(0);

	hsh = &par->HshTab[ HshIdx ];
	hsh->typ1 = typ1;
	hsh->typ2 = typ2;
	mul = pow(typ1->NmbLin, 1./3.) - 1;
	hsh->mul[0] = mul;
	hsh->mul[1] = mul * mul;
	hsh->mul[2] = mul * mul * mul;
	hsh->buc = calloc(typ1->NmbLin, sizeof(BucSct));

	printf("hash mul = %d %d %d\n",hsh->mul[2],hsh->mul[1],hsh->mul[0]);
	printf("hash size = %d, adr = %p\n",typ1->NmbLin,hsh->buc);

	for(i=0;i<par->NmbCpu;i++)
	{
		hsh->ovf[i] = calloc(typ1->NmbLin / par->NmbCpu, sizeof(BucSct));
		hsh->NmbOvf[i] = 0;
	}

	return(HshIdx);
}


/*----------------------------------------------------------*/
/* Release a hash table										*/
/*----------------------------------------------------------*/

void FreeHash(int ParIdx, int HshIdx)
{
	int i;
	ParSct *par;

	if( (ParIdx < 1) || (ParIdx > MaxLibPar) || !(par = ParTab[ ParIdx ]))
		return;

	if( (HshIdx < 1) || (HshIdx > MaxHsh) )
		return;

	if(par->HshTab[ HshIdx ].buc)
		free(par->HshTab[ HshIdx ].buc);

	for(i=0;i<par->NmbCpu;i++)
		if(par->HshTab[ HshIdx ].ovf[i])
			free(par->HshTab[ HshIdx ].ovf[i]);

	memset(&par->HshTab[ HshIdx ], 0, sizeof(HshSct));
}


/*----------------------------------------------------------*/
/* Multithread compatible hash insertertion					*/
/*----------------------------------------------------------*/

long long int AddHash(int ParIdx, int PthIdx, int HshIdx, int a, int b, int c, long long int dat)
{
	int i;
	long long int key, idx[3];
	ParSct *par;
	HshSct *hsh;
	BucSct *buc, *ovf;

	if( (ParIdx < 1) || (ParIdx > MaxLibPar) || !(par = ParTab[ ParIdx ]))
		return(0);

	if( (HshIdx < 1) || (HshIdx > MaxHsh) )
		return(0);

	hsh = &par->HshTab[ HshIdx ];

	if(a < b)
	{
		if(b < c)
		{
			idx[0] = a;
			idx[1] = b;
			idx[2] = c;
		}
		else if(a < c)
		{
			idx[0] = a;
			idx[1] = c;
			idx[2] = b;
		}
		else
		{
			idx[0] = c;
			idx[1] = a;
			idx[2] = b;
		}
	}
	else
	{
		if(a < c)
		{
			idx[0] = b;
			idx[1] = a;
			idx[2] = c;
		}
		else if(b < c)
		{
			idx[0] = b;
			idx[1] = c;
			idx[2] = a;
		}
		else
		{
			idx[0] = c;
			idx[1] = b;
			idx[2] = a;
		}
	}

	key = (hsh->mul[2] * idx[2] + hsh->mul[1] * idx[1] + hsh->mul[0] * idx[0]) / hsh->typ2->NmbLin;
	buc = &hsh->buc[ key ];

/*	printf("pth %d, hash %d / %p, tab = %p, ovf = %p, indices = %lld %lld %lld, key %lld, buc = %p\n",\
		PthIdx, HshIdx, hsh, hsh->buc, hsh->ovf[ PthIdx ], idx[2], idx[1], idx[0], key, buc);
*/
	if(!buc->idx[2])
	{
		for(i=0;i<3;i++)
			buc->idx[i] = idx[i];

		buc->dat = dat;
		return(0);
	}
	else
		do
		{
			if( (buc->idx[0] == idx[0]) && (buc->idx[1] == idx[1]) && (buc->idx[2] == idx[2]) )
				return(buc->dat);
			else if(!buc->nex)
			{
				if(hsh->NmbOvf[ PthIdx ] >= hsh->typ1->NmbLin / par->NmbCpu)
				{
					hsh->ovf[ PthIdx ] = calloc(hsh->typ1->NmbLin / par->NmbCpu, sizeof(BucSct));
					hsh->NmbOvf[ PthIdx ] = 0;
					puts("realloc");
				}

				ovf = &hsh->ovf[ PthIdx ][ hsh->NmbOvf[ PthIdx ]++ ];
				ovf->nex = hsh->buc[ key ].nex;
				hsh->buc[ key ].nex = ovf;
				ovf->dat = dat;

				for(i=0;i<3;i++)
					ovf->idx[i] = idx[i];

				return(0);
			}
		}while(buc = buc->nex);

	return(0);
}


#else


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "lplib3.h"

#define MaxLibPar 10
#define MaxPth 128
#define MaxTyp 100

typedef struct
{
	int NmbLin;
}TypSct;

typedef struct ParSct
{
	int NmbTyp;
	float sta[2];
	void (*prc)(int, int, int, void *), *arg;
	TypSct TypTab[ MaxTyp+1 ];
}ParSct;

ParSct *ParTab[ MaxLibPar+1 ];
int IniLibPar = 0;


int InitParallel(int NmbCpu)
{
	int i, ParIdx;
	ParSct *par = NULL;

	if(NmbCpu > MaxPth)
		return(0);

	if(!IniLibPar)
	{
		IniLibPar = 1;

		for(i=1;i<=MaxLibPar;i++)
			ParTab[i] = NULL;
	}

	/* Allocate and build main parallel structure */

	for(ParIdx=1; ParIdx<=MaxLibPar; ParIdx++)
		if(!ParTab[ ParIdx ])
		{
			par = ParTab[ ParIdx ] = calloc(1, sizeof(ParSct));
			break;
		}

	if(!par)
		return(0);

	return(ParIdx);
}


void StopParallel(int ParIdx)
{
	int i;
	ParSct *par;

	if( (ParIdx < 1) || (ParIdx > MaxLibPar) || !(par = ParTab[ ParIdx ]) )
		return;

	for(i=1;i<=MaxTyp;i++)
		FreeType(ParIdx, i);

	free(par);
	ParTab[ ParIdx ] = NULL;
}


int NewType(int ParIdx, int NmbLin)
{
	int i, TypIdx=0;
	TypSct *typ;
	ParSct *par;

	if( (ParIdx < 1) || (ParIdx > MaxLibPar) || !(par = ParTab[ ParIdx ]) )
		return(0);

	if(NmbLin <= 0)
		return(0);

	for(i=1;i<=MaxTyp;i++)
		if(!par->TypTab[i].NmbLin)
		{
			TypIdx = i;
			break;
		}

	if(!TypIdx)
		return(0);

	typ = &par->TypTab[ TypIdx ];
	typ->NmbLin = NmbLin;

	return(TypIdx);
}


void FreeType(int ParIdx, int TypIdx)
{
}


int BeginDependency(int ParIdx, int TypIdx1, int TypIdx2)
{
	return(0);
}


void AddDependency(int ParIdx, int idx1, int idx2)
{
}


void EndDependency(int ParIdx, float DepSta[2])
{
}


float LaunchParallel(int ParIdx, int typ1, int typ2, void *prc, void *PtrArg)
{
	ParSct *par;
	void (*UsrPrc)(int, int, int, void *) = (void (*)(int, int, int, void *))prc;

	if( (ParIdx < 1) || (ParIdx > MaxLibPar) || !(par = ParTab[ ParIdx ]) )
		return(-1.);

	if( (typ1 < 1) || (typ1 > MaxTyp) || (typ2 < 0) || (typ2 > MaxTyp) || (typ1 == typ2) )
		return(-1.);

	UsrPrc(1, par->TypTab[ typ1 ].NmbLin, 0, PtrArg);

	return(1.);
}


int ParallelMemClear(int ParIdx, void *PtrArg, long siz)
{
	memset(PtrArg, 0, siz);
	return(1);
}

void ParallelQsort(int ParIdx, void *base, size_t nel, size_t width, int (*compar)(const void *, const void *))
{
	ParSct *par;

	/* Get and check lib parallel instance */

	if( (ParIdx < 1) || (ParIdx > MaxLibPar) || !(par = ParTab[ ParIdx ]) )
		return;

	qsort(base, nel, width, compar);
}

#endif
