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
/* SUMMARY : ...                                                            */
/* LICENSE : LGPLv3                                                         */
/* ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE           */
/* AUTHORS : Pascal Frey                                                    */
/* E-MAIL  : pascal.frey@sorbonne-universite.fr                             */

#ifdef __cplusplus
extern "C" {
#endif

#include "medit.h"
#include "extern.h"
#include "sproto.h"
#include "eigenv.h"

int EatLine(FILE *in) {
  int k, c;

  k = 0;

  while ((c = fgetc(in)) != EOF) {
    k++;
    if (c == 10 || c == 13) return (1);
  }

  return (0);
}

int EatSpace(FILE *in) {
  int k, c, ret = 0;

  k = 0;

  while (isspace(c = fgetc(in))) {
    k++;
    if (c == 10 || c == 13) ret = 1;
  }

  if (c == EOF) return -1;

  /*  fprintf(stdout,"EatSpace %d %d last char=%d\n",ret,k,c); */
  ungetc(c, in);

  return ret;
}

int bbfile(pMesh mesh) {
  FILE *in;
  pSolution ps = NULL;
  double a, b, c, lambda[3], eigv[3][3], m[6], vp[2][2];
  float dummy;
  int j, k, l, dim, np, nf, i1, i2, i3, err, iord, ret;
  char *ptr, data[128], tmp[128];
  ubyte bigbb;

  /* default */
  strcpy(tmp, mesh->name);
  ptr = (char *)strstr(tmp, ".mesh");
  if (ptr) *ptr = '\0';

  sprintf(data, "%s.bb", tmp);
  in = fopen(data, "r");
  if (!in) {
    sprintf(data, "%s.pbb", tmp);
    in = fopen(data, "r");
  }

  if (!in) {
    sprintf(data, "%s.BB", tmp);
    in = fopen(data, "r");
    if (!in) { /* hack FH pour le mac */
      sprintf(data, "%s.gbb", tmp);
      in = fopen(data, "r");
    }
  }

  if (!in) return (0);

  i1 = i2 = i3 = -1;
  /* read file format */
  err = 0;
  ret = fscanf(in, "%d", &dim);
  if (ret == EOF) printf("fscanf error\n");
  if (EatSpace(in)) err++;

  ret = fscanf(in, "%d", &i1);
  if (ret == EOF) printf("fscanf error\n");
  if (EatSpace(in)) err++;

  ret = fscanf(in, "%d", &i2);
  if (ret == EOF) printf("fscanf error\n");
  if (EatSpace(in)) err++;

  ret = fscanf(in, "%d", &i3);
  if (ret == EOF) printf("fscanf error\n");
  bigbb = (EatSpace(in) == 0); /* not nl after the 4 integer => BB */

  if (!quiet) {
    if (bigbb)
      fprintf(stdout, "  Reading BB file %s\n", data);
    else
      fprintf(stdout, "  Reading bb file %s\n", data);
  }
  if (dim < 2 || dim > 3 || err) {
    fprintf(stderr, "  %%%% Wrong file (dim=%d) (err=%d). Ignored\n", dim, err);
    return (0);
  }

  /* read number of field(s) */
  nf = 0;

  if (bigbb) {
    int nfield;

    /* get only 1st field */
    nfield = i1;
    mesh->nfield = i2;
    if (nfield > 1) {
      nf += i3;

      for (k = 1; k < nfield - 1; k++) {
        ret = fscanf(in, "%d", &np);
        if (ret == EOF) printf("fscanf error\n");
        nf += np;
      }

      ret = fscanf(in, "%d", &np);
      if (ret == EOF) printf("fscanf error\n");
    } else {
      np = i3;
    }

    /* read file type */
    ret = fscanf(in, "%d", &mesh->typage);
    if (ret == EOF) printf("fscanf error\n");
    printf(" np= %d, type= %d\n", np, mesh->typage);
  } else {
    mesh->nfield = i1;
    np = i2;
    mesh->typage = i3;
  }

  if (mesh->typage == 2) {
    if (np < mesh->np) {
      fprintf(stderr, "  %%%% Wrong solution number (%d , %d). Ignored\n", np, mesh->np);
      fclose(in);
      return (0);
    }

    mesh->nbb = mesh->np;
  } else if (mesh->typage == 1) {
    if (np < mesh->ne) {
      fprintf(stderr, "  %%%% Wrong solution number (%d , %d). Ignored\n", np, mesh->ne);
      fclose(in);
      return (0);
    }

    mesh->nbb = mesh->ne;
  } else {
    fprintf(stderr, "  %%%% Wrong typage (%d). Ignored\n", mesh->typage);
    fclose(in);
    return (0);
  }

  /* read solutions */
  mesh->bbmin = 1.e10;
  mesh->bbmax = -1.e10;

  /* allocate memory */
  if (!zaldy2(mesh)) {
    mesh->nbb = 0;
    fclose(in);
    return (0);
  }

  /* scalar field */
  if (mesh->nfield == 1) {
    if (ddebug) printf("   scalar (isotropic) field\n");

    for (k = 1; k <= mesh->nbb; k++) {
      ps = &mesh->sol[k];
      ps->bb = 0.0;
      if (fscanf(in, "%127s", data) != 1) continue;

      if ((ptr = strpbrk(data, "dD"))) *ptr = 'E';

      ret = sscanf(data, "%f", &ps->bb);
      if (ret == EOF) printf("sscanf error\n");
      if (ps->bb < mesh->bbmin) mesh->bbmin = ps->bb;

      if (ps->bb > mesh->bbmax) mesh->bbmax = ps->bb;

      for (j = 1; j <= nf; j++) {
        ret = fscanf(in, "%f", &dummy);
        if (ret == EOF) printf("sscanf error\n");
      }
    }
  }

  /* vector field */
  else if (mesh->nfield == mesh->dim) {
    if (ddebug) fprintf(stdout, "   vector field \n");

    for (k = 1; k <= mesh->nbb; k++) {
      ps = &mesh->sol[k];
      ps->bb = 0.0;

      for (l = 0; l < mesh->dim; l++) {
        if (fscanf(in, "%127s", data) != 1) continue;

        if ((ptr = strpbrk(data, "dD"))) *ptr = 'E';

        ret = sscanf(data, "%f", &ps->m[l]);
        if (ret == EOF) printf("sscanf error\n");
        ps->bb += ps->m[l] * ps->m[l];
      }

      ps->bb = sqrt(ps->bb);
      if (ps->bb < mesh->bbmin) mesh->bbmin = ps->bb;

      if (ps->bb > mesh->bbmax) mesh->bbmax = ps->bb;

      for (j = 1; j < nf; j++) {
        ret = fscanf(in, "%f", &dummy);
        if (ret == EOF) printf("fscanf error\n");
      }
    }

    fclose(in);
    return (0);
  } else if (dim == 2 && mesh->nfield == 3) {
    if (ddebug) fprintf(stdout, "   2D metric field\n");

    for (k = 1; k <= mesh->np; k++) {
      ps = &mesh->sol[k];
      ret = fscanf(in, "%lf %lf %lf", &a, &b, &c);
      if (ret == EOF) printf("fscanf error\n");
      ps->m[0] = a;
      ps->m[1] = b;
      ps->m[2] = c;
      m[0] = a;
      m[1] = b;
      m[2] = c;
      eigen2(m, lambda, vp);
      ps->bb = min(lambda[0], lambda[1]);
      if (ps->bb < mesh->bbmin) mesh->bbmin = ps->bb;

      if (ps->bb > mesh->bbmax) mesh->bbmax = ps->bb;

      for (j = 1; j < nf; j++) {
        ret = fscanf(in, "%f", &dummy);
        if (ret == EOF) printf("fscanf error\n");
      }
    }
  } else if (dim == 3 && mesh->nfield == 6) {
    if (ddebug) fprintf(stdout, "   3D metric field\n");

    for (k = 1; k <= mesh->np; k++) {
      ps = &mesh->sol[k];
      ps->bb = 0.0f;

      for (l = 0; l < 6; l++) {
        if (fscanf(in, "%127s", data) != 1) {
          continue;
          m[l] = 0;
        }

        if ((ptr = strpbrk(data, "dD"))) *ptr = 'E';

        ret = sscanf(data, "%f", &dummy);
        if (ret == EOF) printf("sscanf error\n");
        m[l] = dummy;
      }

      ps->m[0] = m[0];
      ps->m[1] = m[1];
      ps->m[2] = m[3];
      ps->m[3] = m[2];
      ps->m[4] = m[4];
      ps->m[5] = m[5];

      m[2] = ps->m[2];
      m[3] = ps->m[3];
      iord = eigenv(1, m, lambda, eigv);
      if (iord) {
        ps->bb = lambda[0];
        ps->bb = max(ps->bb, lambda[1]);
        ps->bb = max(ps->bb, lambda[2]);
        if (ps->bb < mesh->bbmin) mesh->bbmin = ps->bb;

        if (ps->bb > mesh->bbmax) mesh->bbmax = ps->bb;
      } else {
        fprintf(stdout, "  ## Eigenvalue problem.\n");
      }

      for (j = 1; j < nf; j++) {
        ret = fscanf(in, "%f", &dummy);
        if (ret == EOF) printf("fscanf error\n");
      }
    }
  } else {
    fprintf(stderr, " %%%% Solution not suitable. Ignored\n");
    mesh->nbb = 0;
  }

  fclose(in);
  return (np);
}

#ifdef __cplusplus
}
#endif
