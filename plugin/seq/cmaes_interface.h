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
// SUMMARY : CMA-ES for non-linear function minimization
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : Nikolaus Hansen
// E-MAIL  : ...

#ifndef _CMAES_INTERFACE_H_
#define _CMAES_INTERFACE_H_

#include "cmaes.h"

/* --------------------------------------------------------- */
/* ------------------ Interface ---------------------------- */
/* --------------------------------------------------------- */

/* --- initialization, constructors, destructors --- */
double *cmaes_init(cmaes_t *, int dimension, double *xstart, double *stddev, long seed, int lambda,
                   const char *input_parameter_filename);
void cmaes_resume_distribution(cmaes_t *evo_ptr, char *filename);
void cmaes_exit(cmaes_t *);

/* --- core functions --- */
double *const *cmaes_SamplePopulation(cmaes_t *);
double *cmaes_UpdateDistribution(cmaes_t *, const double *rgFitnessValues);
const char *cmaes_TestForTermination(cmaes_t *);

/* --- additional functions --- */
double *const *cmaes_ReSampleSingle(cmaes_t *t, int index);
double const *cmaes_ReSampleSingle_old(cmaes_t *, double *rgx);
double *cmaes_SampleSingleInto(cmaes_t *t, double *rgx);
void cmaes_UpdateEigensystem(cmaes_t *, int flgforce);

/* --- getter functions --- */
double cmaes_Get(cmaes_t *, char const *keyword);
const double *cmaes_GetPtr(cmaes_t *, char const *keyword); /* e.g. "xbestever" */
double *cmaes_GetNew(cmaes_t *t, char const *keyword);
double *cmaes_GetInto(cmaes_t *t, char const *keyword, double *mem);

/* --- online control and output --- */
void cmaes_ReadSignals(cmaes_t *, char const *filename);
void cmaes_WriteToFile(cmaes_t *, const char *szKeyWord, const char *output_filename);
char *cmaes_SayHello(cmaes_t *);
/* --- misc --- */
double *cmaes_NewDouble(int n);
void cmaes_FATAL(char const *s1, char const *s2, char const *s3, char const *s4);

#endif
