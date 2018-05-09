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

#include "medit.h"
#include "extern.h"
#include "sproto.h"

int zaldy2 (pMesh mesh) {
	pSolution ps;
	int k, nbf;

	/* memory alloc. */
	mesh->sol = (pSolution)M_calloc(mesh->nbb + 1, sizeof(struct solu), "zaldy2");
	assert(mesh->sol);

	if (mesh->nfield == 1) return (1);

	if (mesh->nfield == mesh->dim)
		nbf = mesh->dim;/* vector field */
	else
		nbf = mesh->dim * (mesh->dim + 1) / 2;	/* d*d matrix */

	for (k = 1; k <= mesh->nbb; k++) {
		ps = &mesh->sol[k];
		ps->m = (float *)malloc(nbf * sizeof(float));
	}

	return (1);
}
