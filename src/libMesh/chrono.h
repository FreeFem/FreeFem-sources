/*
 * This file is part of FreeFem++.
 *
 * FreeFem++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * FreeFem++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifdef __cplusplus
extern "C" {
#endif

#include <time.h>

#ifndef  ON
#define  RESET 0
#define  ON 1
#define  OFF 2
#endif

#define  TIMEMAX 16
#define  MAXCLK (1073741823. / (double)CLOCKS_PER_SEC)

typedef struct mytime {
	double ctim, dtim;
	time_t ptim;
	short call;
} mytime;

/* prototypes */
void chrono (int cmode, mytime *ptt);
double gttime (mytime t);
void tminit (mytime *t, int maxtim);

#ifdef __cplusplus
}
#endif
