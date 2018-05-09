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

#include <stdlib.h>
#include <assert.h>

/* prototype (re)definitions */
void*M_malloc (size_t size, char *call);
void*M_calloc (size_t nelem, size_t elsize, char *call);
void*M_realloc (void *ptr, size_t size, char *call);
void M_free (void *ptr);

/* ptototypes : tools */
int M_memLeak ();
void M_memDump ();
size_t M_memSize ();

#ifdef __cplusplus
}
#endif
