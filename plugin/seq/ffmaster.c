/*
 * Example of coupling c program and freemfem++ script
 * with mmap and semaphore
 *
 * the c code is    :   ffmaster.c
 * the ff++ code is : ffslave.edp
 * and here FreeFem++ is a slave process
 * the compile step is
 *
 * cc -c libff-mmap-semaphore.c
 * cc ffmaster.c -o ffmaster  libff-mmap-semaphore.o -g
 #build the freefem++ plugin
 * ff-c++ -auto ff-mmap-semaphore.cpp
 # launch
 # ./ffmaster
 #
 #
 # F. Hecht Feb. 2018   Frederic.Hecht@upmc.fr
 */

/*
 *
 * This file is part of Freefem++
 *
 * Freefem++ is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2.1 of the License, or
 * (at your option) any later version.
 *
 * Freefem++  is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Freefem++; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "libff-mmap-semaphore.h"
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
ff_Psem sem_ff, sem_c;	// the semaphore for mutex
/*
 * Psemaphore smff("ff-slave");
 * Psemaphore smc("ff-master");
 * Pmmap sharedata("shared-data");
 *
 */
int main (int argc, const char **argv) {
	int debug = 0;
	ff_Pmmap shd;
	double cff, rff;
	long status;
	int i, ret;

	if (argc > 1) {debug = atoi(argv[1]);}

	ff_mmap_sem_verb = debug;

	sem_ff = ffsem_malloc();
	sem_c = ffsem_malloc();
	shd = ffmmap_malloc();

	ffsem_init(sem_ff, "ff-slave1", 1);
	ffsem_init(sem_c, "ff-master1", 1);
	ffmmap_init(shd, "shared-data", 1024);

	status = 1;
	ffmmap_write(shd, &status, sizeof(status), 8);
	ffmmap_msync(shd, 0, 32);

	char ff[1024];
	sprintf(ff, "../../src/nw/FreeFem++ ../../examples/plugin/ffslave.edp -nw -ns -v %d&", debug);
	ret = system(ff);	// Lauch FF++ in batch no graphique
	if (ret == -1) printf("system function error\n");
	if (debug) {printf(" cc: before wait\n");}

	if (debug) {printf(" cc: before wait 0 ff\n");}

	ffsem_wait(sem_ff);

	for (i = 0; i < 10; ++i) {
		printf(" iter : %d \n", i);
		cff = 10 + i;
		ffmmap_write(shd, &cff, sizeof(cff), 0);
		ffsem_post(sem_c);

		if (debug) {printf(" cc: before wait 2\n");}

		ffsem_wait(sem_ff);
		ffmmap_read(shd, &rff, sizeof(rff), 16);
		printf(" iter = %d rff= %f\n", i, rff);
	}

	status = 0;	// Fin
	ffmmap_write(shd, &status, sizeof(status), 8);
	ffsem_post(sem_c);
	printf("Fin Master \n");
	ffsem_wait(sem_ff);
	ffsem_del(sem_ff);
	ffsem_del(sem_c);
	ffmmap_del(shd);
	return 0;
}
