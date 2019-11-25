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

#define FLOAT_MAX (float)1.e20

int loadGIS(pMesh mesh) {
  pQuad pq;
  pPoint ppt;
  pSolution ps;
  FILE *fp;
  double xxm, yym, ggx, ggy, hhz;
  float *te, cx, cy, cz, gu, hu, xmi, ymi;
  int i, j, k, sx, sy, ret, bitsize, ref;
  char *ptr, c, buf[256], data[256];
  ubyte ityp;

  /* default */
  strcpy(data, mesh->name);
  ptr = strstr(data, ".gis");
  if (!ptr) strcat(data, ".gis");

  fp = fopen(data, "rb");
  if (!fp) return (0);

  fprintf(stdout, "  Reading %s\n", data);

  if (!fgets(buf, sizeof(buf), fp)) {
    fprintf(stderr, "  ## Invalid header.\n");
    fclose(fp);
    return (0);
  }

  /* remove leading spaces */
  /* 13: carriage return */
  while (buf[0] == ' ') {    //  && ((buf[0] != 13) && (buf[0] != '\n'))) {
    memmove(&buf[0], &buf[1], strlen(buf));
  }

  /* check header file */
  if (buf[0] != 'G') {
    fprintf(stderr, "  ## Invalid format.\n");
    fclose(fp);
    return (0);
  }

  if (buf[1] == '1') {
    ityp = 1;
  } else if (buf[1] == '2') {
    ityp = 2;
  } else {
    fprintf(stderr, "  ## Invalid format ('G?' expected).\n");
    fclose(fp);
    return (0);
  }

  /* check and strip comments */
  do {
    ret = fscanf(fp, "%255s", buf);
    if (ret == EOF) break;

    if (buf[0] == '#') do {
        c = getc(fp);
      } while (c != '\n');

    else
      break;
  } while (1);

  /* header */
  ret = sscanf(buf, "%d", &sx);
  ret += fscanf(fp, "%d", &sy);
  ret += fscanf(fp, "%f %f %f", &cx, &cy, &cz);
  ret += fscanf(fp, "%f", &gu);
  ret += fscanf(fp, "%f", &hu);
  ret += fscanf(fp, "%f %f", &xmi, &ymi);
  if (ret != 9) {
    fprintf(stderr, "  ## Error loading terrain.\n");
    free(mesh);
    fclose(fp);
    return (0);
  }

  bitsize = mesh->np * sizeof(float);
  if (ddebug) {
    fprintf(stdout, "   header: sx  %7d   sy %7d\n", sx, sy);
    fprintf(stdout, "           cx  %7.2f   cy %7.2f  cz %7.2f\n", cx, cy, cz);
    fprintf(stdout, "        units: %7.2f  %7.2f\n", gu, hu);
    fprintf(stdout, "          min: %7.3e  %7.3e\n", xmi, ymi);
    fprintf(stdout, "   terrain size: %dx%d  %ld bytes\n", sx, sy, (long)bitsize);
  }

  mesh->np = sx * sy;
  mesh->nt = 0;
  mesh->nq = (sx - 1) * (sy - 1);
  mesh->ne = mesh->nq;

  /* memory allocation for mesh */
  if (!zaldy1(mesh)) {
    fclose(fp);
    return (0);
  }

  /* read data */
  if (ityp == 1) {
    xxm = xmi / cx;
    yym = ymi / cy;
    ggx = gu * cx;
    ggy = gu * cy;
    hhz = hu * cz;

    for (j = 1; j <= sy; j++) {
      for (i = 1; i <= sx; i++) {
        k = (j - 1) * sx + i;
        ppt = &mesh->point[k];
        ret = fscanf(fp, "%lf", &ppt->c[2]);
        if (ret == EOF) printf("fgets error\n");
        ppt->c[0] = (float)(ggx * (xxm + i - 1));
        ppt->c[1] = (float)(ggy * (yym + j - 1));
        ppt->c[2] = (float)(hhz * ppt->c[2]);
      }
    }
  } else {
    int ni;

    te = (float *)malloc(sx * sizeof(float));
    if (!te) exit(1);

    ni = 0;
    xxm = xmi / cx;
    yym = ymi / cy;
    ggx = gu / cx;
    ggy = gu / cy;
    hhz = hu / cz;

    for (j = 1; j <= sy; j++) {
      ret = fread(te, sizeof(int), sx, fp);
      if (ret != sx) {
        fprintf(stderr, "  ## Error loading terrain.\n");
        free(mesh->point);
        free(mesh);
        free(te);
        fclose(fp);
        return (0);
      }

      for (i = 1; i <= sx; i++) {
        k = (j - 1) * sx + i;
        ppt = &mesh->point[k];
        ppt->c[2] = te[ni++];
        ppt->c[0] = ggx * (xxm + i - 1);
        ppt->c[1] = ggy * (yym + j - 1);
        ppt->c[2] = hhz * ppt->c[2];
      }
    }

    free(te);
  }

  /* construct topology */
  mesh->dim = 3;

  for (j = 1; j < sy; j++) {
    for (i = 1; i < sx; i++) {
      k = (j - 1) * (sx - 1) + i;
      pq = &mesh->quad[k];
      pq->v[0] = (j - 1) * sx + i;
      pq->v[1] = pq->v[0] + 1;
      pq->v[2] = j * sx + i + 1;
      pq->v[3] = pq->v[2] - 1;
      pq->ref = 0;
    }
  }

  /* read references, if any */
  if (ityp == 1) {
    int pos = ftell(fp);
    ret = fscanf(fp, "%d", &i);
    if (ret != EOF) {
      fseek(fp, pos, SEEK_SET);

      for (j = 1; j < sy; j++) {
        for (i = 1; i < sx; i++) {
          k = (j - 1) * (sx - 1) + i;
          pq = &mesh->quad[k];
          ret = fscanf(fp, "%d", &ref);
          if (ret == EOF) printf("fgets error\n");
          pq->ref = ref;
        }
      }
    }
  }

  /* solution = elevations */
  mesh->nbb = mesh->np;
  mesh->bbmin = FLOAT_MAX;
  mesh->bbmax = -FLOAT_MAX;
  mesh->typage = 2;
  mesh->nfield = 1;

  /* allocate memory */
  if (!zaldy2(mesh)) {
    mesh->nbb = 0;
    fclose(fp);
    return (1);
  }

  for (j = 1; j <= mesh->np; j++) {
    ps = &mesh->sol[j];
    ppt = &mesh->point[j];
    ps->bb = ppt->c[2];
    if (ps->bb < mesh->bbmin) mesh->bbmin = ps->bb;

    if (ps->bb > mesh->bbmax) mesh->bbmax = ps->bb;
  }

  fclose(fp);
  return (1);
}

#ifdef __cplusplus
}
#endif
