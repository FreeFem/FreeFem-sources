// Emacs will be in -*- Mode: c++ -*-
//
// ********** DO NOT REMOVE THIS BANNER **********
//
// SUMMARY: Language for a Finite Element Method
// RELEASE: 2.0     
// USAGE  : You may copy freely these files and use it for    
//          teaching or research. These or part of these may   
//          not be sold or used for a commercial purpose with- 
//          out our consent : fax (33)1 44 27 44 11        
//
// AUTHORS:  D. Bernardi, Y. Darmaillac F. Hecht,    
//           P. Parole O. Pironneau C. Prud'homme
// ORG    :          
// E-MAIL :   pironneau@ann.jussieu.fr     
//
// ORIG-DATE:     June-94
// LAST-MOD:     16-Jan-96 at 16:10:19 by Prud'homme Christophe
//
// DESCRIPTION:  
// DESCRIP-END.
//

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
namespace Fem2D {
#define i_nint(x) ((long) x >= 0 ? x + 0.5 : x - 0.5)
#define amin(a, b) ((a) < (b) ? (a) : (b))
#define amax(a, b) ((a) < (b) ? (b) : (a))
#define aabs(a) ((a) < 0 ? -(a) : (a))

static long     nothing = -1073741824L;



/*****************************************************************************
I. Fabrique une triangulation (de nbsmax sommets maxi.) a partir de:
--------------------------------------------------------------------
    -cr: tableau de taille 2*nbs (abscisse & ordonnees des points en entree)
    -arete: tableau de taille 2*nba (depart & arrivee des aretes)
      (decrire le contour dans le sens direct: domaine a
            gauche des aretes exterieures)
    -h: tableau de taille nbs (precision relative aux points en entree)


II. Alloue et affecte une structure t de type triangulation:
------------------------------------------------------------
     -t.np: nombre de sommets en sortie
     -t.nt: nombre de triangles en sortie
     -t.rp[i].x:} abscisse
      t.rp[i].y:} & ordonnee du ieme point du maillage
     -t.tr[i][j]: numero du jeme point du ieme triangle du maillage


III. Renvoie le type d'erreur:
------------------------------
     -  0 : pas d'erreur
     -(-1): pas assez de memoire
     - >0 : erreur mshptg_
  mshptg erreurs
1: mshptg := 'le nb de point est < 3 or > que le nb max de point';
2: mshptg := 'il y des points confondus ';
3: mshptg := 'tout les points sont alignªs ';
7: mshptg := 'bug dans l''algorithme pour retrouver la frontiüre (Internal Bug mshfrt)';
8: mshptg := 'Une arçte a forcer traverse trop de triangles (Internal Bug : mshfr1)';
9: mshptg := 'la frontiüre est croisªe ';
10: mshptg := 'il y a un point qui est sur une arçte frontiüre';
11: mshptg := 'Un sous domaine est reference par aucun triangle ';
20: mshptg := 'mshopt: 3 points confondus (Internal Bug) ';
21: mshptg := 'la pile de mshopt est trop petit (Internal Bug)';



    Remarque: penser a liberer les pointeurs t.rp et t.tr!!!
    --------- un tableau de taille N est indexe de 0 a N-1
        les numeros affectes aux points et aux triangles commencent a 0
******************************************************************************/
/*
long 
MakeTriangulation (triangulation * t, long nbs, long nbsmax, long nba,
		   float *crbdy, float *hbdy, long *arete, int *ngbdy, long *sd, long nbsd, int* flag, int fflag)
{
  int             i, j;
  long            err = 0, nbt, nbsold = nbs;
  long           *c = NULL;
  long           *tri = NULL;
  long           *nu = NULL;
  long           *reft = NULL;
  int            *ngg = NULL;
  float          *cr = NULL;
  float          *h = NULL;

  nbt = 2 * nbsmax;
  nu = new long[6*nbt];
  c = new long[2*nbsmax];
  ngg = new int[nbsmax];
  tri = new long[(4 * nbsmax + 2 * nbsd)];
  reft = new long[nbt];
  cr = new float[(2 * nbsmax + 2)];
  h = new float[nbsmax];

  for (i = 0; i < 2 * nba; i++)
    arete[i]++;
  for (i = 0; i < nbs; i++)
     {
       ngg[i] = ngbdy[i];
       cr[2 * i] = crbdy[2 * i];
       cr[2 * i + 1] = crbdy[2 * i + 1];
       h[i] = hbdy[i];
     }
  for (i = nbs; i < nbsmax; i++)
    ngg[i] = 0;

  mshptg_ (cr, h, c, nu, &nbs, nbsmax, tri, arete, nba, (long *) sd, nbsd, reft, &nbt, .25, .75, &err);
  for (i = 0; i < 2 * nba; i++)
    arete[i]--;
  if (err)
    goto L1;
  if (*flag)
     {
       delete [] t->rp;t->rp = NULL;
       delete [] t->tr;t->tr = NULL;
       delete [] t->ng;t->ng = NULL;
       delete [] t->ngt;t->ngt = NULL;
     }
  t->rp =new rpoint[nbs];
  t->tr = new triangle[nbt];
  t->ng = new int[nbs];
  t->ngt = new int[nbt];

  *flag = 1;
  t->np = nbs;
  t->nt = nbt;
  for (i = 0; i < nbt; i++)
     {
       for (j = 0; j < 3; j++)
	 t->tr[i][j] = nu[3 * i + j] - 1;
       t->ngt[i] = reft[i] - 1;
     }
  for (i = 0; i < nbs; i++)
     {
       t->rp[i].x = cr[2 * i];
       t->rp[i].y = cr[2 * i + 1];
       if (i < nbsold)
	 t->ng[i] = ngg[i];
       else
	 t->ng[i] = 0;
     }
  renum (t);
  if(!fflag) removeBdyT (t);
L1:
  delete [] nu;nu = NULL;
  delete [] cr;cr = NULL;
  delete [] c;c = NULL;
  delete [] tri;tri = NULL;
  delete [] reft;reft = NULL;
  delete [] ngg;ngg  = NULL;
  delete [] h;h = NULL;
  return err;
}



int 
verifietr (float *cr, int ns)
{
  int             i;
  float           xh, diameter = 1.F;

  if (ns == 0)
    return -1;
  if (ns > 1)
     {
       diameter = 0.F;
       for (i = 0; i < ns; i++)
	 diameter = amax (diameter, aabs (cr[2 * i] - cr[0]) + aabs (cr[2 * i + 1] - cr[1]));
     }
  else
    diameter = 0.001F;
  for (i = 0; i < ns; i++)
     {
       xh = aabs (cr[2 * i] - cr[2 * ns]) + aabs (cr[2 * i + 1] - cr[2 * ns + 1]);
       if (xh < 1e-5 * diameter)
	 return i;
     }
  return -1;
}

*/

int             mshrgl_ (float *c, long *nrfs, long *nbs, long *nu, long *w1,
			 long *w, float omega, long itermx, float eps);
int             mshopt_ (long *c, long *nu, long *t, long a, long *err);
void            mshvoi_ (long *nu, long *w1, long *w, long *nbt, long *nbs);
int             msha1p_ (long *t, long *s, long *c, long *nu, long *reft, long *tete, long *nbt,
			 long *err);
int             mshtri_ (float *cr, long *c, long *nbs, long *tri, long *nu, float *trfri, long *err);
int             mshcxi_ (long *c, long *nu, long *tri, long *nbs, long *tete, long *err);
int             mshfrt_ (long *c, long *nu, long *nbs, long *arete, long nba, long *sd,
			 long nbsd, long *reft, long *w, long *err);
int             mshgpt_ (long *c, float *cr, long *nu, float *h, long *reft, long *nbs,
    long nbsmx, long *nbt, float coef, float puis, float *trfri, long *err);
long            mshlcl_ (long *c, long *nu, long *tete, long *s);
int             mshtr1_ (long *criter, long *record, long *n);
int             mshcvx_ (long direct, long *c, long *nu, long *pfold, long *err);
int             mshfr1_ (long *c, long *nu, long *it1, long *ita, long *is1, long *s2, long *err);
int             mshfr2_ (long *c, long *nu, long *lst, long *nbac, long *t, long *ta,
			 long *ss1, long *ss2, long *err);

int 
mshptg_ (float *cr, float *h, long *c, long *nu, long *nbs, long nbsmx, long *tri,
	 long *arete, long nba, long *sd,
	 long nbsd, long *reft, long *nbt, float coef, float puis, long *err)
{
  /* System generated locals */
  long            i_1;

  /* Local variables */
  static long     tete, i, j, k, t;
  static float    trfri[4];
  static long     nbsgrn, nbtgrn;

/* ----------------------------------------------------------------------- */
/*      but:  construire une triangulation a partir d'un ensemble de */
/*             points et d'un maillage frontalier */
/* ----------------------------------------------------------------------- */
/*     entre : */
/*     ------- */
/*           cr(2,nbsmx)  tableau des coordonnees des nbs points donnes */


/*           h (nbsmx)    tableau du h local voulu autour de chaque point */
/*                          donnes */
/*           nbs          nombre de points donnes */
/*           nbsmx        nombre de points maximal a cree */
/*                        si nbs  = nbsmx ont ne cree pas de points */
/*                        si nbs  < nbsmx => erreur */
/*           arete(2,nba) tableau des aretes du maillage a forcer */
/*                          exemple :la frontiere */
/*           nba          le nombre d'aretes du maillage */
/*           sd(2,nbsd)   tableau definisant les nbsd  sous domaine */
/*                          (reference des triangles gerener) */
/*                          abs(sd(1,i)) =  numero d'une l'arete */
/*                          si sd(1,i) est positive alors le sous domaine */
/*                          est a gauche de l'arete sinon il est a droite */
/*                          sd(2,i) donne le numero du sous doimaine */

/*           puis         coefficent de generation des points */
/*                        .1  => on propage plus loin les rafinement */
/*                               donnes par h */
/*                        .25 => valeur conseillee */
/*           coef         coefficent sur le test arret */
/*                          le valeur conseillee est .75 */
/*                          remarque le nombre de point genere est en */
/*                          O(coef**2) */

/*        tableaux de travail: */
/*        -------------------- */
/*           c(2,nbsmx)    tableau d'entiers (copie de coordonnees) */
/*           tri(ltri)     tableau d'entiers */
/*        out : */
/*        ----- */
/*         nbs         nombre de points   donnes + generes */
/*         nbt         nombre de triangles generes */
/*         cr(1:2,nbs) coordonnees des sommets donnes + generes */
/*         nu(1:3,nbt) sommets des triangles (tableau des connections) */
/*                       telle que les sommets tourne dans le sens direct */
/*         reft(1:nbt) numero de sous domaine de chaque triangle */
/*         err    si err = 0 alors pas de probleme */
/*                sinon nbt = 0 et pas de triangulation */
/*     dimension des tableaux */
/*     ---------------------- */
/*     definition des parameters */
/*     nbtmx = 2*(nbs-1) ,  ltri = max(4*nbs+2*nbsd,nba) */
/*     long : nu(6*nbtmx) , reft(nbtmx) , c(2*nbsmx) , tri(ltri) */
/*     long : arete(2,nba), sd(2,nbsd) */
/*     float    : cr(2*nbsmx) , h(nbsmx) */

/* ---------------------------------------------------------------------- */
/*     programmeur F.Hecht, France */
/*           version 1.0  mars 1986 */
/* ----------------------------------------------------------------------- */

  /* Parameter adjustments */
  --reft;
  sd -= 3;
  arete -= 3;
  --tri;
  --nu;
  c -= 3;
  --h;
  cr -= 3;

  /* Function Body */
  *err = 0;
  *nbt = 0;
  if (*nbs < 3 || nbsmx < *nbs)
     {
       *err = 1;
       return 0;
     }
/* ------------------------- */
/* preparation des donnees */
/* ------------------------- */
  mshtri_ (&cr[3], &c[3], nbs, &tri[1], &tri[*nbs + 1], trfri, err);
  if (*err != 0)
     {
       return 0;
     }
/* -------------------------------- */
/* maillage de l enveloppe convexe */
/* -------------------------------- */
  mshcxi_ (&c[3], &nu[1], &tri[1], nbs, &tete, err);
/* ----------------------------------------------------------------------- */
/*     definition de tableau nu(1:6,2*nbs-2) */
/* ----------------------------------------------------------------------- */
/*     nu(*,ie) definit soit un element ,soit un sommet frontiere */
/*     si nu(5:6,ie) = (0,0) alors ie est un sommet frontiere */
/*     avec nu(1,ie) = numero du sommet */
/*          nu(2,ie) = 8*t + a */
/*                     ou t est le numero du triangle ayant l'arete */
/*                     frontiere (a) dont le premier sommet est nu(1,ie) */
/*          nu(3,ie) = pointeur dans nu sur sommet frontiere precedent */
/*          nu(4,ie) = pointeur dans nu sur sommet frontiere suivant */
/*     sinon ie est un element : */
/*          nu(1:3,ie) numero des 3 sommets du triangle ie tournant dans */
/*                     le sens direct */
/*          nu(4:6,ie) = (d4,d5,d6) donnee des 3 aretes ai */
/*           ai est forme des sommets nu(i-3,ie),nu(mod(i,3)+1,ie) */
/*           si di < 0 alors arete i est frontiere et -di est pointeur */
/*             sur 1er sommet frontiere de i */
/*           sinon arete est interne et di = 8*ta + ata */
/*              ou ta est le numero du triangle adjacent a l'arete */
/*              et ata est le numero de l'arete dans ta */
/* ------------------------------------------------------------------------ */
  if (*err != 0)
     {
       return 0;
     }

  i_1 = *nbs;
  for (i = 1; i <= i_1; ++i)
     {
       tri[i] = 0;
     }
  i = tete;
L20:
  j = nu[(i - 1) * 6 + 4];
  tri[nu[(i - 1) * 6 + 1]] = nu[(j - 1) * 6 + 1];
  i = j;
  if (i != tete)
     {
       goto L20;
     }
/* ----------------------------- */
/* traitement frontiere */
/* ----------------------------- */
  k = 0;
  mshfrt_ (&c[3], &nu[1], nbs, &arete[3], nba, &sd[3], nbsd, &reft[1], &tri[
								   1], err);
  if (*err != 0)
     {
       return 0;
     }
/* ------------------------------------------------------------------- */
/*       on a modifie nu les sommets frontiere n'ont plus de sens */
/*       ainsi que les pointeurs sur ces elements */
/* ------------------------------------------------------------------- */
  nbsgrn = *nbs;
  mshgpt_ (&c[3], &cr[3], &nu[1], &h[1], &reft[1], &nbsgrn, nbsmx, &nbtgrn,
	   coef, puis, trfri, err);
  if (*err != 0)
     {
       return 0;
     }
/*     construction du tableau nu(1:3,1:nbt) */
/* ------------------------------------------ */
  *nbt = 0;
  k = 0;
  j = 1;
  i_1 = nbtgrn;
  for (t = 1; t <= i_1; ++t)
     {
       if (nu[j + 5] != 0)
	  {
	    ++(*nbt);
	    reft[*nbt] = reft[t];
	    for (i = 0; i <= 2; ++i)
	       {
		 ++k;
		 nu[k] = nu[j + i];
	       }
	  }
       j += 6;
     }
/*     dans nu il y a (s1(t),s2(t),s3(t),t=1,nbt) */
/*     ou s1 s2 s3 sont les 3 sommets de t */
/* ------------------------------------------------ */
  i_1 = *nbs;
  for (i = 1; i <= i_1; ++i)
     {
       tri[i] = 1;
     }
  i_1 = nbsgrn;
  for (i = *nbs + 1; i <= i_1; ++i)
     {
       tri[i] = 0;
     }
  mshvoi_ (&nu[1], &tri[nbsgrn + 1], &nu[*nbt * 3 + 1], nbt, &nbsgrn);
  mshrgl_ (&cr[3], &tri[1], &nbsgrn, &nu[1], &tri[nbsgrn + 1], &nu[*nbt * 3
						      + 1], 1.4F, 20L, .005F);
  *nbs = nbsgrn;
  return 1;
}				/* mshptg_ */

/* ********************************************************************** */
void 
mshvoi_ (long *nu, long *w1, long *w, long *nbt, long *nbs)
{
  /* System generated locals */
  long            i_1;

  /* Local variables */
  static long     i, is;


/* RECHERCHE DU VOISINAGE */
/* ----------------------- */
  /* Parameter adjustments */
  --w;
  --nu;

  /* Function Body */
  i_1 = *nbs;
  for (i = 1; i <= i_1; ++i)
     {
       w1[i] = 0;
     }
  i_1 = *nbt * 3;
  for (i = 1; i <= i_1; ++i)
     {
       ++w1[nu[i]];
     }
  w1[0] = 0;
  i_1 = *nbs;
  for (i = 1; i <= i_1; ++i)
     {
       w1[i] = w1[i - 1] + w1[i];
     }
  i_1 = *nbt * 3;
  for (i = 1; i <= i_1; ++i)
     {
       is = nu[i] - 1;
       ++w1[is];
       w[w1[is]] = i;
     }
  for (i = *nbs; i >= 1; --i)
     {
       w1[i] = w1[i - 1];
     }
  w1[0] = 0;
}				/* mshvoi_ */

/* ********************************************************************** */
int 
mshrgl_ (float *c, long *nrfs, long *nbs, long *nu, long *w1,
	 long *w, float omega, long itermx, float eps)
{
  /* System generated locals */
  long            i_1, i_2, i_3;
  float           r_1, r_2;

  /* Local variables */
  static float    depx, depy;
  static long     iter;
  static float    xmin, ymin, xmax, ymax;
  static long     i, k, i1, i2, ic;
  static float    bx, by, dx;
  static long     is;
  static float    err;


/* REGULARISATION PAR MOYENNE BARYCENTRIQUE */
/* ----------------------------------------- */
  /* Parameter adjustments */
  --w;
  --nu;
  --nrfs;
  c -= 3;

  /* Function Body */
  xmin = c[3];
  ymin = c[4];
  xmax = c[3];
  ymax = c[4];
  i_1 = *nbs;
  for (ic = 1; ic <= i_1; ++ic)
     {
/* Computing MIN */
       r_1 = c[(ic << 1) + 1];
       xmin = amin (r_1, xmin);
/* Computing MIN */
       r_1 = c[(ic << 1) + 2];
       ymin = amin (r_1, ymin);
/* Computing MAX */
       r_1 = c[(ic << 1) + 1];
       xmax = amax (r_1, xmax);
/* Computing MAX */
       r_1 = c[(ic << 1) + 2];
       ymax = amax (r_1, ymax);
     }
/* Computing MAX */
  r_1 = xmax - xmin, r_2 = ymax - ymin;
  dx = amax (r_1, r_2);
  i_1 = itermx;
  for (iter = 1; iter <= i_1; ++iter)
     {
       err = (float) 0.;
       i2 = w1[0];
       i_2 = *nbs;
       for (is = 1; is <= i_2; ++is)
	  {
	    i1 = i2 + 1;
	    i2 = w1[is];
	    if (i2 >= i1 && nrfs[is] == 0)
	       {
		 bx = (float) 0.;
		 by = (float) 0.;
		 i_3 = i2;
		 for (i = i1; i <= i_3; ++i)
		    {
		      if (w[i] % 3 == 0)
			 {
			   k = w[i] - 2;
			 }
		      else
			 {
			   k = w[i] + 1;
			 }
		      bx += c[(nu[k] << 1) + 1];
		      by += c[(nu[k] << 1) + 2];
		    }
		 bx /= i2 - i1 + 1;
		 by /= i2 - i1 + 1;
		 depx = omega * (c[(is << 1) + 1] - bx);
		 depy = omega * (c[(is << 1) + 2] - by);
		 c[(is << 1) + 1] -= depx;
		 c[(is << 1) + 2] -= depy;
/* Computing MAX */
		 r_1 = err, r_2 = aabs (depx), r_1 = amax (r_1, r_2), r_2 = aabs (depy);
		 err = amax (r_1, r_2);
	       }
	  }
/* -------------------------------- */
       if (err <= eps * dx)
	  {
	    return 0;
	  }
     }
  return 1;
}				/* mshrgl_ */

/* ********************************************************************** */
int 
mshgpt_ (long *c, float *cr, long *nu, float *h, long *reft, long *nbs,
     long nbsmx, long *nbt, float coef, float puis, float *trfri, long *err)
{
  /* System generated locals */
  long            i_1;
  float           r_1, r_2, r_3;
  double          d_1, d_2, d_3, d_4, d_5, d_6, d_7, d_8;

  /* Local variables */
  static float    aire;
  static long     tete;
  static long     t;
  static float    x, y;
  static long     itera;
  static float    h1, h2, h3;
  static long     s1, s2, s3;
  static float    hs;
  static long     ix, iy, nbsold;
  static float    det, pui;

  /* Parameter adjustments */
  --trfri;
  --reft;
  --h;
  nu -= 7;
  cr -= 3;
  c -= 3;

  /* Function Body */
  pui = puis;
  *nbt = (*nbs << 1) - 2;
  if (*nbs >= nbsmx)
     {
       return 0;
     }
  tete = 0;
/*     initialisation de la liste des triangles libre */
  i_1 = *nbt;
  for (t = 1; t <= i_1; ++t)
     {
       if (nu[t * 6 + 6] == 0)
	  {
	    nu[t * 6 + 1] = tete;
	    tete = t;
	  }
     }
  itera = 0;
L20:
  ++itera;
  nbsold = *nbs;
  i_1 = *nbt;
  for (t = 1; t <= i_1; ++t)
     {
       if (nu[t * 6 + 6] != 0)
	  {
	    s1 = nu[t * 6 + 1];
	    s2 = nu[t * 6 + 2];
	    s3 = nu[t * 6 + 3];
/*        calcul de 2 fois l'aire du triangle */
	    det = (cr[(s2 << 1) + 1] - cr[(s1 << 1) + 1]) * (cr[(s3 << 1) + 2]
		    - cr[(s1 << 1) + 2]) - (cr[(s2 << 1) + 2] - cr[(s1 << 1)
			    + 2]) * (cr[(s3 << 1) + 1] - cr[(s1 << 1) + 1]);
	    aire = det * coef;
	    if (puis > (float) 0.)
	       {
		 d_2 = (double) h[s1];
		 d_3 = (double) pui;
		 d_4 = (double) h[s2];
		 d_5 = (double) pui;
		 d_6 = (double) h[s3];
		 d_7 = (double) pui;
		 d_1 = (pow (d_2, d_3) + pow (d_4, d_5) + pow (d_6, d_7)) / (float) 3.;
		 d_8 = (double) (1.F / pui);
		 hs = (float)pow (d_1, d_8);
	       }
	    else if (puis > (float) -1.)
	       {
		 d_1 = (double) (h[s1] * h[s2] * h[s3]);
		 hs = (float)pow (d_1, 1. / 3);
	       }
	    else if (puis > (float) -2.)
	       {
		 hs = h[s1] * (float) 3. *h[s2] * h[s3] / (h[s1] * h[s2] + h[
					       s1] * h[s3] + h[s2] * h[s3]);
	       }
	    else
	       {
/* Computing 2nd power */
		 r_1 = (float)sqrt (h[s1] * h[s2]);
/* Computing 2nd power */
		 r_2 = h[s1] * h[s3];
/* Computing 2nd power */
		 r_3 = h[s2] * h[s3];
		 hs = (float)sqrt (3.0) * (h[s1] * h[s2] * h[s3]
			 / (r_1 * r_1) + r_2 * r_2 + r_3 * r_3);
	       }
	    if (aire > hs * hs)
	       {
		 h1 = (float) 1.;
		 h2 = (float) 1.;
		 h3 = (float) 1.;
		 x = (cr[(s1 << 1) + 1] * h1 + cr[(s2 << 1) + 1] * h2 + cr[(s3
					  << 1) + 1] * h3) / (h1 + h2 + h3);
		 y = (cr[(s1 << 1) + 2] * h1 + cr[(s2 << 1) + 2] * h2 + cr[(s3
					  << 1) + 2] * h3) / (h1 + h2 + h3);
		 r_1 = trfri[1] * (x - trfri[2]);
		 ix = (long) i_nint (r_1);
		 r_1 = trfri[1] * (y - trfri[3]) - trfri[4];
		 iy = (long) i_nint (r_1);
		 if ((c[(s2 << 1) + 1] - ix) * (c[(s3 << 1) + 2] - iy) - (c[(
		       s2 << 1) + 2] - iy) * (c[(s3 << 1) + 1] - ix) <= 0 ||
		     (ix - c[(s1 << 1) + 1]) * (c[(s3 << 1) + 2] - c[(s1 <<
			 1) + 2]) - (iy - c[(s1 << 1) + 2]) * (c[(s3 << 1) +
		      1] - c[(s1 << 1) + 1]) <= 0 || (c[(s2 << 1) + 1] - c[(
			s1 << 1) + 1]) * (iy - c[(s1 << 1) + 2]) - (c[(s2 <<
		       1) + 2] - c[(s1 << 1) + 2]) * (ix - c[(s1 << 1) + 1])
		     <= 0)
		    {
		    }
		 else
		    {
		      if (*nbs >= nbsmx)
			 {
			   return 0;
			 }
		      ++(*nbs);
		      c[(*nbs << 1) + 1] = ix;
		      c[(*nbs << 1) + 2] = iy;
		      cr[(*nbs << 1) + 1] = ix / trfri[1] + trfri[2];
		      cr[(*nbs << 1) + 2] = (iy + trfri[4]) / trfri[1] + trfri[
									 3];
		      h[*nbs] = hs;
		      msha1p_ (&t, nbs, &c[3], &nu[7], &reft[1], &tete, nbt, err);
		      if (*err != 0)
			 {
			   return 0;
			 }
		    }
	       }
	  }
     }
  if (nbsold != *nbs)
     {
       goto L20;
     }
  return 1;
}				/* mshgpt_ */

/* ********************************************************************** */
int 
msha1p_ (long *t, long *s, long *c, long *nu, long *reft, long *tete, long *nbt,
	 long *err)
{
  static long     t1, t2, t3, ia2, ia3;
  static long     ta2, ta3;

//    extern int mshopt_();
  static long     tta;

  /* Parameter adjustments */
  --reft;
  nu -= 7;
  c -= 3;

  /* Function Body */
  t1 = *t;
  if (*tete == 0)
     {
       ++(*nbt);
       t2 = *nbt;
     }
  else
     {
       t2 = *tete;
       *tete = nu[*tete * 6 + 1];
     }
  if (*tete == 0)
     {
       ++(*nbt);
       t3 = *nbt;
     }
  else
     {
       t3 = *tete;
       *tete = nu[*tete * 6 + 1];
     }
  nu[t2 * 6 + 1] = *s;
  nu[t2 * 6 + 2] = nu[*t * 6 + 2];
  nu[t2 * 6 + 3] = nu[*t * 6 + 3];
  nu[t2 * 6 + 4] = (t1 << 3) + 5;
  nu[t2 * 6 + 5] = nu[*t * 6 + 5];
  nu[t2 * 6 + 6] = (t3 << 3) + 5;
  nu[t3 * 6 + 1] = nu[*t * 6 + 1];
  nu[t3 * 6 + 2] = *s;
  nu[t3 * 6 + 3] = nu[*t * 6 + 3];
  nu[t3 * 6 + 4] = (t1 << 3) + 6;
  nu[t3 * 6 + 5] = (t2 << 3) + 6;
  nu[t3 * 6 + 6] = nu[*t * 6 + 6];
  tta = nu[*t * 6 + 5];
  if (tta > 0)
     {
       ta2 = tta / 8;
       ia2 = tta - (ta2 << 3);
       nu[ia2 + ta2 * 6] = (t2 << 3) + 5;
     }
  tta = nu[*t * 6 + 6];
  if (tta > 0)
     {
       ta3 = tta / 8;
       ia3 = tta - (ta3 << 3);
       nu[ia3 + ta3 * 6] = (t3 << 3) + 6;
     }
  nu[t1 * 6 + 3] = *s;
  nu[t1 * 6 + 5] = (t2 << 3) + 4;
  nu[t1 * 6 + 6] = (t3 << 3) + 4;
  reft[t2] = reft[*t];
  reft[t3] = reft[*t];
  mshopt_ (&c[3], &nu[7], &t1, 4L, err);
  if (*err != 0)
     {
       return 0;
     }
  mshopt_ (&c[3], &nu[7], &t2, 5L, err);
  if (*err != 0)
     {
       return 0;
     }
  mshopt_ (&c[3], &nu[7], &t3, 6L, err);
  if (*err != 0)
     {
       return 0;
     }
  return 1;
}				/* msha1p_ */

/* ********************************************************************** */
long 
mshlcl_ (long *c, long *nu, long *tete, long *s)
{
  /* System generated locals */
  long            ret_val;

  /* Local variables */
  static long     init;
  static long     x, y, pt, det, ppt;

  /* Parameter adjustments */
  nu -= 7;
  c -= 3;

  /* Function Body */
  x = c[(*s << 1) + 1];
  y = c[(*s << 1) + 2];
  init = 1;
  pt = *tete;
L10:
  ppt = pt;
  pt = nu[pt * 6 + 4];
  if (pt != *tete)
     {
       det = x * c[(nu[pt * 6 + 1] << 1) + 2] - y * c[(nu[pt * 6 + 1] << 1)
						      + 1];
       if (det < 0)
	  {
	    init = 0;
	    goto L10;
	  }
       else if (init && det == 0)
	  {
	    goto L10;
	  }
     }
  ret_val = ppt;
  return ret_val;
}				/* mshlcl_ */

/* ********************************************************************** */
int 
mshtri_ (float *cr, long *c, long *nbs, long *tri, long *nu, float *trfri, long *err)
{
  /* System generated locals */
  long            i_1, i_2, i_3;
  float           r_1;

  /* Local variables */
  static long     ierr, trik;
  static float    xmin, ymin, xmax, ymax;
  static long     i, j, k, ic, jc;

//    extern int mshtr1_();
  static long     ip, xx;
  static float    aa1, aa2;
  static long     iii, det, tri3;

  /* Parameter adjustments */
  --trfri;
  --nu;
  --tri;
  c -= 3;
  cr -= 3;

  /* Function Body */
  ierr = 0;
  iii = 1;
  xmin = cr[3];
  ymin = cr[4];
  xmax = cr[3];
  ymax = cr[4];
  i_1 = *nbs;
  for (ic = 1; ic <= i_1; ++ic)
     {
/* Computing MIN */
       r_1 = cr[(ic << 1) + 1];
       xmin = amin (r_1, xmin);
/* Computing MIN */
       r_1 = cr[(ic << 1) + 2];
       ymin = amin (r_1, ymin);
/* Computing MAX */
       r_1 = cr[(ic << 1) + 1];
       xmax = amax (r_1, xmax);
/* Computing MAX */
       r_1 = cr[(ic << 1) + 2];
       ymax = amax (r_1, ymax);
       tri[ic] = ic;
       if (cr[(ic << 1) + 1] < cr[(iii << 1) + 1])
	  {
	    iii = ic;
	  }
     }
  aa1 = (float) 32767. / (xmax - xmin);
  aa2 = (float) 32767. / (ymax - ymin);
  aa1 = amin (aa1, aa2);
  aa2 = aa1 * (cr[(iii << 1) + 2] - ymin);
  trfri[1] = aa1;
  trfri[2] = cr[(iii << 1) + 1];
  trfri[3] = ymin;
  trfri[4] = aa2;
  i_1 = *nbs;
  for (ic = 1; ic <= i_1; ++ic)
     {
       r_1 = aa1 * (cr[(ic << 1) + 1] - cr[(iii << 1) + 1]);
       c[(ic << 1) + 1] = (long) i_nint (r_1);
       r_1 = aa1 * (cr[(ic << 1) + 2] - ymin) - aa2;
       c[(ic << 1) + 2] = (long) i_nint (r_1);
/* Computing 2nd power */
       i_2 = c[(ic << 1) + 1];
/* Computing 2nd power */
       i_3 = c[(ic << 1) + 2];
       nu[ic] = i_2 * i_2 + i_3 * i_3;
     }
/* ---------------------------------------------------------- */
  mshtr1_ (&nu[1], &tri[1], nbs);
  ip = 1;
  xx = nu[ip];
  i_1 = *nbs;
  for (jc = 1; jc <= i_1; ++jc)
     {
       if (nu[jc] > xx)
	  {
	    i_2 = jc - ip;
	    mshtr1_ (&nu[ip], &tri[ip], &i_2);
	    i_2 = jc - 2;
	    for (i = ip; i <= i_2; ++i)
	       {
		 if (nu[i] == nu[i + 1])
		    {
		      ++ierr;
		    }
	       }
	    xx = nu[jc];
	    ip = jc;
	  }
       ic = tri[jc];
       nu[jc] = c[(ic << 1) + 2];
     }
  i_1 = *nbs - ip;
  mshtr1_ (&nu[ip], &tri[ip], &i_1);
  i_1 = jc - 2;
  for (i = ip; i <= i_1; ++i)
     {
       if (nu[i] == nu[i + 1])
	  {
	    ++ierr;
	  }
     }
  if (ierr != 0)
     {
       *err = 2;
       return 0;
     }
  k = 2;
L50:
  if (k <= *nbs)
     {
       ++k;
       det = c[(tri[2] << 1) + 1] * c[(tri[k] << 1) + 2] - c[(tri[2] << 1) +
						  2] * c[(tri[k] << 1) + 1];
       if (det == 0)
	  {
	    goto L50;
	  }
     }
  else
     {
       *err = 3;
       return 0;
     }
/*     k est le premier point non aligne */
  trik = tri[k];
  for (j = k - 1; j >= 3; --j)
     {
       tri[j + 1] = tri[j];
     }
  tri[3] = trik;
  if (det < 0)
     {
/*       on inverse les  points 2 3 tries */
       tri3 = tri[3];
       tri[3] = tri[2];
       tri[2] = tri3;
     }
  return 1;
}				/* mshtri_ */

/* ********************************************************************** */
int 
mshtr1_ (long *criter, long *record, long *n)
{
  /* System generated locals */
  long            i_1;

  /* Local variables */
  static long     crit, i, j, l, r, rec;

/*     trie selon les valeurs de criter croissantes */
/*     record suit le reordonnancement */


  /* Parameter adjustments */
  --record;
  --criter;

  /* Function Body */
  if (*n <= 1)
     {
       return 0;
     }
  l = *n / 2 + 1;
  r = *n;
L2:
  if (l <= 1)
     {
       goto L20;
     }
  --l;
  rec = record[l];
  crit = criter[l];
  goto L3;
L20:
  rec = record[r];
  crit = criter[r];
  record[r] = record[1];
  criter[r] = criter[1];
  --r;
  if (r == 1)
     {
       goto L999;
     }
L3:
  j = l;
L4:
  i = j;
  j <<= 1;
  if ((i_1 = j - r) < 0)
     {
       goto L5;
     }
  else if (i_1 == 0)
     {
       goto L6;
     }
  else
     {
       goto L8;
     }
L5:
  if (criter[j] < criter[j + 1])
     {
       ++j;
     }
L6:
  if (crit >= criter[j])
     {
       goto L8;
     }
  record[i] = record[j];
  criter[i] = criter[j];
  goto L4;
L8:
  record[i] = rec;
  criter[i] = crit;
  goto L2;
L999:
  record[1] = rec;
  criter[1] = crit;
  return 0;
}				/* mshtr1_ */

/* ********************************************************************** */
int 
mshcvx_ (long direct, long *c, long *nu, long *pfold, long *err)
{
  static long     t, a4, a5, i1, i2, i3, i4, i5, i6, s1, s2, s3, t4, t5,
                  pf, pp, ps;

//    extern int mshopt_();
  static long     tt4, tt5, det, ppf, psf;

  /* Parameter adjustments */
  nu -= 7;
  c -= 3;

  /* Function Body */
  if (direct)
     {
       pp = 3;
       ps = 4;
       i1 = 1;
       i2 = 3;
       i3 = 2;
       i4 = 6;
       i5 = 5;
       i6 = 4;
     }
  else
     {
       pp = 4;
       ps = 3;
       i1 = 1;
       i2 = 2;
       i3 = 3;
       i4 = 4;
       i5 = 5;
       i6 = 6;
     }
L10:
  ppf = *pfold;
  pf = nu[ps + *pfold * 6];
  psf = nu[ps + pf * 6];
  s1 = nu[ppf * 6 + 1];
  s2 = nu[pf * 6 + 1];
  s3 = nu[psf * 6 + 1];
  det = (c[(s2 << 1) + 1] - c[(s1 << 1) + 1]) * (c[(s3 << 1) + 2] - c[(s1 <<
	     1) + 2]) - (c[(s2 << 1) + 2] - c[(s1 << 1) + 2]) * (c[(s3 << 1)
						   + 1] - c[(s1 << 1) + 1]);
  if (!(direct) && det > 0 || direct && det < 0)
     {
/*       on ajoute un triangle t et on detruit une arete */
/*       ----------------------------------------------- */
       if (direct)
	  {
	    tt4 = nu[ppf * 6 + 2];
	    tt5 = nu[pf * 6 + 2];
	  }
       else
	  {
	    tt4 = nu[pf * 6 + 2];
	    tt5 = nu[psf * 6 + 2];
	  }
       t4 = tt4 / 8;
       t5 = tt5 / 8;
       a4 = tt4 - (t4 << 3);
       a5 = tt5 - (t5 << 3);
/*       destruction de l'arete frontiere en pf */
/*       -------------------------------------- */
       nu[ps + ppf * 6] = psf;
       nu[pp + psf * 6] = ppf;
/*       on remplace l'arete frontiere par l'element genere */
/*       --------------------------------------------------- */
       t = pf;
/*       update de l'arete non detruite */
/*       ------------------------------ */
       if (direct)
	  {
	    nu[ppf * 6 + 2] = (t << 3) + i6;
	  }
       else
	  {
	    nu[psf * 6 + 2] = (t << 3) + i6;
	  }
/*       on cree l'element */
/*       ----------------- */
       nu[i1 + t * 6] = s1;
       nu[i2 + t * 6] = s2;
       nu[i3 + t * 6] = s3;
       nu[i4 + t * 6] = (t4 << 3) + a4;
       nu[i5 + t * 6] = (t5 << 3) + a5;
       if (direct)
	  {
	    nu[i6 + t * 6] = -ppf;
	  }
       else
	  {
	    nu[i6 + t * 6] = -psf;
	  }
       nu[a4 + t4 * 6] = (t << 3) + i4;
       nu[a5 + t5 * 6] = (t << 3) + i5;
       mshopt_ (&c[3], &nu[7], &t5, a5, err);
       if (*err != 0)
	  {
	    return 0;
	  }
       goto L10;
     }
  return 1;
}				/* mshcvx_ */

/* ********************************************************************** */
int 
mshcxi_ (long *c, long *nu, long *tri, long *nbs, long *tete, long *err)
{
  /* System generated locals */
  long            i_1;

  /* Local variables */
  static long     sfree, ttaf, i, j, s, t, pf;

//    extern long mshlcl_();
  //    extern int mshcvx_(), mshopt_();
  static long     iaf, taf, npf, ppf, psf;

/*     initialisation de la sfree liste dans nu */
  /* Parameter adjustments */
  --tri;
  nu -= 7;
  c -= 3;

  /* Function Body */
  i_1 = *nbs + *nbs - 2;
  for (i = 1; i <= i_1; ++i)
     {
       nu[i * 6 + 1] = i + 1;
       for (j = 2; j <= 6; ++j)
	  {
	    nu[j + i * 6] = 0;
	  }
     }
  nu[(*nbs + *nbs - 2) * 6 + 1] = 0;
  sfree = 1;
/*     initialisation du premier triangle */
  t = sfree;
  sfree = nu[sfree * 6 + 1];
/*     initialisation de la liste frontiere */
  *tete = sfree;
  pf = sfree;
  for (i = 1; i <= 3; ++i)
     {
       nu[i + t * 6] = tri[i];
       nu[i + 3 + t * 6] = -pf;
       ppf = pf;
       sfree = nu[pf * 6 + 1];
       pf = sfree;
       if (i == 3)
	  {
	    pf = *tete;
	  }
       nu[ppf * 6 + 1] = tri[i];
       nu[ppf * 6 + 2] = i + 3 + (t << 3);
       nu[ppf * 6 + 4] = pf;
       nu[pf * 6 + 3] = ppf;
     }
  i_1 = *nbs;
  for (i = 4; i <= i_1; ++i)
     {
       s = tri[i];
       pf = mshlcl_ (&c[3], &nu[7], tete, &s);
/*      creation d'un nouveau triangle et modification de la
   frontiere */
/*      
   --------------------------------------------------------------
 */
       t = sfree;
       sfree = nu[sfree * 6 + 1];
       npf = sfree;
       sfree = nu[sfree * 6 + 1];
       ppf = nu[pf * 6 + 3];
       psf = nu[pf * 6 + 4];
       ttaf = nu[pf * 6 + 2];
       taf = ttaf / 8;
       iaf = ttaf - (taf << 3);

/*                  npf */
/*               1  x s               --- */
/*                 / \                --- */
/*              4 /   \ 6        ---  vide --- */
/*               /  t  \              --- */
/*            2 /   5   \ 3           --- */
/* ------ --<---x---------x---------x- frontiere--<--- */
/*          psf \  iaf  /  pf         --- */
/*               \ taf /         --- omega --- */
/*                \   /               --- */
/*                 \ /                --- */
/*                  x                 --- */
/*                                    --- */
/*     generation  de l'element t */
       nu[t * 6 + 1] = s;
       nu[t * 6 + 2] = nu[psf * 6 + 1];
       nu[t * 6 + 3] = nu[pf * 6 + 1];
       nu[t * 6 + 4] = -npf;
       nu[t * 6 + 5] = (taf << 3) + iaf;
       nu[t * 6 + 6] = -pf;
       nu[iaf + taf * 6] = (t << 3) + 5;
/*      update de la liste frontiere */
       nu[npf * 6 + 4] = psf;
       nu[pf * 6 + 4] = npf;
       nu[npf * 6 + 3] = pf;
       nu[psf * 6 + 3] = npf;
       nu[npf * 6 + 1] = s;
       nu[npf * 6 + 2] = (t << 3) + 4;
       nu[pf * 6 + 2] = (t << 3) + 6;
       mshopt_ (&c[3], &nu[7], &t, 5L, err);
       if (*err != 0)
	  {
	    return 0;
	  }
       mshcvx_ (1, &c[3], &nu[7], &npf, err);
       if (*err != 0)
	  {
	    return 0;
	  }
       mshcvx_ (0, &c[3], &nu[7], &npf, err);
       if (*err != 0)
	  {
	    return 0;
	  }
     }
  return 1;
}				/* mshcxi_ */

/* ********************************************************************** */
int 
mshopt_ (long *c, long *nu, long *t, long a, long *err)
{
  /* Initialized data */

  static long     mod3[3] =
  {2, 3, 1};

  /* System generated locals */
  long            i_1;
  double          d_1;

  /* Local variables */
  static long     pile[4096] /* was [2][256] */ ;
  static float    reel1, reel2;
  static double   reel8;
  static long     i, a1, a2, s1, t1, t2, s2, s3, s4, aa, i11, i12, i13,
                  i21, i22, i23, tt;
  static long     tt1, sgn, cos1, cos2, sin1, sin2;

  /* Parameter adjustments */
  nu -= 7;
  c -= 3;

  /* Function Body */
  i = 1;
  pile[(i << 1) - 2] = *t;
  pile[(i << 1) - 1] = a;
L10:
  if (i > 0)
     {
       t1 = pile[(i << 1) - 2];
       a1 = pile[(i << 1) - 1];
       --i;
       if (t1 <= 0)
	  {
	    goto L10;
	  }
       tt1 = nu[a1 + t1 * 6];
       if (tt1 <= 0)
	  {
	    goto L10;
	  }
       t2 = tt1 / 8;
       a2 = tt1 - (t2 << 3);
       i11 = a1 - 3;
       i12 = mod3[i11 - 1];
       i13 = mod3[i12 - 1];
       i21 = a2 - 3;
       i22 = mod3[i21 - 1];
       i23 = mod3[i22 - 1];
       s1 = nu[i13 + t1 * 6];
       s2 = nu[i11 + t1 * 6];
       s3 = nu[i12 + t1 * 6];
       s4 = nu[i23 + t2 * 6];
       sin1 = (c[(s3 << 1) + 2] - c[(s1 << 1) + 2]) * (c[(s2 << 1) + 1] - c[(
	       s1 << 1) + 1]) - (c[(s3 << 1) + 1] - c[(s1 << 1) + 1]) * (c[(
					  s2 << 1) + 2] - c[(s1 << 1) + 2]);
       cos1 = (c[(s3 << 1) + 1] - c[(s1 << 1) + 1]) * (c[(s3 << 1) + 1] - c[(
	       s2 << 1) + 1]) + (c[(s3 << 1) + 2] - c[(s1 << 1) + 2]) * (c[(
					  s3 << 1) + 2] - c[(s2 << 1) + 2]);
       if (sin1 == 0 && cos1 == 0)
	  {
	    *err = 20;
	    return 0;
	  }
/*       b est la cotangente de angle (s1,s3,s2) */
       sin2 = (c[(s4 << 1) + 1] - c[(s1 << 1) + 1]) * (c[(s2 << 1) + 2] - c[(
	       s1 << 1) + 2]) - (c[(s4 << 1) + 2] - c[(s1 << 1) + 2]) * (c[(
					  s2 << 1) + 1] - c[(s1 << 1) + 1]);
       cos2 = (c[(s4 << 1) + 1] - c[(s2 << 1) + 1]) * (c[(s4 << 1) + 1] - c[(
	       s1 << 1) + 1]) + (c[(s4 << 1) + 2] - c[(s2 << 1) + 2]) * (c[(
					  s4 << 1) + 2] - c[(s1 << 1) + 2]);
       reel1 = (float) cos2 *(float) sin1;
       reel2 = (float) cos1 *(float) sin2;

       if (aabs (reel1) + aabs (reel2) >= (float) 1073741824.)
	  {
	    reel8 = (double) cos2 *(double) sin1 + (double) cos1 *(double) sin2;

/* Computing MIN */
	    d_1 = amax (reel8, -1.);
	    reel8 = amin (d_1, 1.);
	    sgn = (long) reel8;
	  }
       else
	  {
	    sgn = cos2 * sin1 + cos1 * sin2;
	  }
/* Computing MIN */
       i_1 = amax (sgn, -1);
       if (amin (i_1, 1) * sin1 >= 0)
	  {
	    goto L10;
	  }
/*       on inverse le quadrilatere */
/*       update des sommets */
/* ------------------------- */
       nu[i12 + t1 * 6] = s4;
       nu[i22 + t2 * 6] = s1;
/*       update des aretes a1,a2 */
/* ------------------------------- */
       tt1 = nu[i22 + 3 + t2 * 6];
       nu[a1 + t1 * 6] = tt1;
       if (tt1 > 0)
	  {
	    tt = tt1 / 8;
	    aa = tt1 - (tt << 3);
	    nu[aa + tt * 6] = a1 + (t1 << 3);
	  }
       else if (tt1 != nothing)
	  {
	    nu[-tt1 * 6 + 2] = a1 + (t1 << 3);
	  }
       tt1 = nu[i12 + 3 + t1 * 6];
       nu[a2 + t2 * 6] = tt1;
       if (tt1 > 0)
	  {
	    tt = tt1 / 8;
	    aa = tt1 - (tt << 3);
	    nu[aa + tt * 6] = a2 + (t2 << 3);
	  }
       else if (tt1 != nothing)
	  {
	    nu[-tt1 * 6 + 2] = a2 + (t2 << 3);
	  }
       nu[i12 + 3 + t1 * 6] = i22 + 3 + (t2 << 3);
       nu[i22 + 3 + t2 * 6] = i12 + 3 + (t1 << 3);
       if (i + 4 > 256)
	  {
	    *err = 21;
	    return 0;
	  }
       ++i;
       pile[(i << 1) - 2] = t1;
       pile[(i << 1) - 1] = a1;
       ++i;
       pile[(i << 1) - 2] = t2;
       pile[(i << 1) - 1] = a2;
       ++i;
       pile[(i << 1) - 2] = t1;
       pile[(i << 1) - 1] = i13 + 3;
       ++i;
       pile[(i << 1) - 2] = t2;
       pile[(i << 1) - 1] = i23 + 3;
       goto L10;
     }
  return 1;
}				/* mshopt_ */

/* ********************************************************************** */
int 
mshfrt_ (long *c, long *nu, long *nbs, long *arete, long nba, long *sd,
	 long nbsd, long *reft, long *w, long *err)
{
  /* Initialized data */

  static long     p3[3] =
  {2, 3, 1};

  /* System generated locals */
  long            i_1, i_2, i_3;

  /* Local variables */
  static long     nbac, ifrt, a, i, t, itera, s1, s2;

//    extern int mshfr1_();
  static long     ie, ap, ta, is, nbacpp;
  static long     is1, ss1, s2t, s3t, isd, jsd, nbt, det2, det3, err1;

  /* Parameter adjustments */
  --w;
  --reft;
  sd -= 3;
  arete -= 3;
  nu -= 7;
  c -= 3;

  /* Function Body */
  if (nba == 0)
     {
       return 0;
     }
  ifrt = 0;
  nbt = *nbs + *nbs - 2;
  i_1 = *nbs;
  for (i = 1; i <= i_1; ++i)
     {
       reft[i] = 0;
     }
  i_1 = nba;
  for (i = 1; i <= i_1; ++i)
     {
       reft[arete[(i << 1) + 1]] = nothing;
       reft[arete[(i << 1) + 2]] = nothing;
     }
  nbac = 0;
  i_1 = nba;
  for (a = 1; a <= i_1; ++a)
     {
/* Computing MIN */
       i_2 = arete[(a << 1) + 1], i_3 = arete[(a << 1) + 2];
       s1 = amin (i_2, i_3);
/* Computing MAX */
       i_2 = arete[(a << 1) + 1], i_3 = arete[(a << 1) + 2];
       s2 = amax (i_2, i_3);
       if (s1 == s2)
	  {
	    ++nbac;
	  }
       else
	  {
	    i = reft[s1];
	  L25:
	    if (i != nothing)
	       {
/* Computing MAX */
		 i_2 = arete[(i << 1) + 1], i_3 = arete[(i << 1) + 2];
		 if (s2 == amax (i_2, i_3))
		    {
		      ++nbac;
		    }
		 else
		    {
		      i = w[i];
		      goto L25;
		    }
	       }
	    else
	       {
		 w[a] = reft[s1];
		 reft[s1] = a;
	       }
	  }
     }
  nbacpp = 1;
  itera = 0;
  err1 = 0;
L50:
  ++itera;
  if (err1 != 0)
     {
       *err = err1;
       return 0;
     }
  if (nbac < nba)
     {
       if (nbacpp == 0)
	  {
	    i_1 = *nbs;
	    for (i = 1; i <= i_1; ++i)
	       {
		 a = reft[i];
	       L60:
		 if (a > 0)
		    {
		      s1 = arete[(i << 1) + 1];
		      s2 = arete[(i << 1) + 2];
		      a = w[a];
		      goto L60;
		    }
	       }
	    *err = 7;
	    return 0;
	  }
/* --------------------------------------------------------------------- */
/*     on s'occupe des aretes a forcer */
/* --------------------------------------------------------------------- */
       nbacpp = 0;
       i_1 = nbt;
       for (ie = 1; ie <= i_1; ++ie)
	  {
	    if (nu[ie * 6 + 5] != 0)
	       {
		 for (is = 1; is <= 3; ++is)
		    {
		      s1 = nu[is + ie * 6];
		      s2t = nu[p3[is - 1] + ie * 6];
		      ss1 = amin (s1, s2t);
		      ap = 0;
		      a = reft[ss1];
		    L80:
		      if (a > 0)
			 {
/* Computing MAX */
			   i_2 = arete[(a << 1) + 1], i_3 = arete[(a << 1) + 2];
			   s2 = amax (i_2, i_3);
			   t = ie;
			   ta = 0;
			   if (s2 == amax (s1, s2t))
			      {
				if (nu[is + 3 + ie * 6] > 0)
				   {
				     ta = nu[is + 3 + ie * 6] / 8;
				     i = nu[is + 3 + ie * 6] - (ta << 3);
				     nu[i + ta * 6] = nothing;
				   }
				nu[is + 3 + ie * 6] = nothing;
				goto L100;
			      }
			   ap = a;
			   a = w[a];
			   goto L80;
			 }
		      if (itera == 1)
			 {
			   goto L110;
			 }
		      ss1 = s1;
		      ap = 0;
		      a = reft[ss1];
		    L90:
		      if (a > 0)
			 {
/* Computing MAX */
			   i_2 = arete[(a << 1) + 1], i_3 = arete[(a << 1) + 2];
			   s2 = amax (i_2, i_3);
			   t = ie;
			   ta = 0;
/*             recherche si l' element coupe l''arete a */
			   is1 = is;
			   s3t = nu[p3[p3[is - 1] - 1] + t * 6];
			   det2 = (c[(s2t << 1) + 1] - c[(s1 << 1) + 1]) * (c[(
			      s2 << 1) + 2] - c[(s1 << 1) + 2]) - (c[(s2t <<
				1) + 2] - c[(s1 << 1) + 2]) * (c[(s2 << 1) +
						     1] - c[(s1 << 1) + 1]);
			   det3 = (c[(s3t << 1) + 1] - c[(s1 << 1) + 1]) * (c[(
			      s2 << 1) + 2] - c[(s1 << 1) + 2]) - (c[(s3t <<
				1) + 2] - c[(s1 << 1) + 2]) * (c[(s2 << 1) +
						     1] - c[(s1 << 1) + 1]);
			   if (det2 > 0 && det3 < 0)
			      {
				mshfr1_ (&c[3], &nu[7], &t, &ta, &is1, &s2, err);
				if (*err != 0)
				   {
				     return 0;
				   }
				goto L100;
			      }
			   else if (det2 == 0 && reft[s2t] == 0)
			      {
				err1 = 10;
			      }
			   else if (det3 == 0 && reft[s3t] == 0)
			      {
				err1 = 10;
			      }
			   ap = a;
			   a = w[a];
			   goto L90;
			 }
		      goto L110;
		    L100:
		      ++nbacpp;
		      if (ap == 0)
			 {
			   reft[ss1] = w[a];
			 }
		      else
			 {
			   w[ap] = w[a];
			 }
		      if (nbac + nbacpp == nba)
			 {
			   goto L130;
			 }
		    L110:
		      ;
		    }
	       }
	  }
       nbac += nbacpp;
       goto L50;
     }
L130:
/* ----------------------------------------------------------------------- */
/*     prise en compte des sous domaines */
/* ----------------------------------------------------------------------- */
/*  add FH   si pas de nbsd ---  jan 2004 */
    if (nbsd == 0) {
	 long i__1 = nbt,nbssd,exter,i__,headt,nst,j;
	for (t = 1; t <= i__1; ++t) {
	    reft[t] = -1073741824;
	}
	nbssd = 0;
/*         if(impre.gt.1) print *,'nbsd .eq. 0 =>  recherche de ssd' */
	i__1 = nbt;
	for (t = 1; t <= i__1; ++t) {
	    if (nu[t * 6 + 1] > 0 && nu[t * 6 + 5] != 0) {
		++nbssd;
		exter = 0;
		i__ = 2;
		w[i__ - 1] = t;
		w[i__] = 3;
		reft[t] = 0;
		headt = t;
		nu[t * 6 + 1] = -nu[t * 6 + 1];
		nst = 1;
/*              print *,' ssd ',nbssd */
L131:
		if (i__ > 0) {
		    ++w[i__];
		    if (w[i__] <= 6) {
			ta = nu[w[i__] + w[i__ - 1] * 6];
			if (ta <= 0) {
			    if (ta != -1073741824) {
				exter = 1;
			    }
			} else if (ta > 0) {
			    ta /= 8;
			    if (nu[ta * 6 + 1] > 0) {
/*                           print *,ta */
				nu[ta * 6 + 1] = -nu[ta * 6 + 1];
				++nst;
				reft[ta] = headt;
				headt = ta;
				w[i__ + 1] = ta;
				w[i__ + 2] = 3;
				i__ += 2;
			    }
			}
		    } else {
			i__ += -2;
		    }
		    goto L131;
		}

		if (exter) {
		    --nbssd;
		    i__ = headt;
L133:
		    if (i__ > 0) {
			j = reft[i__];
/*                     print *,i */
			reft[i__] = 0;
			nu[i__ * 6 + 1] = 0;
			i__ = j;
			goto L133;
		    }
		} else {
		    i__ = headt;
L136:
		    if (i__ > 0) {
			j = reft[i__];
			reft[i__] = nbssd;
			i__ = j;
			goto L136;
		    }
		}
	    }
	}
	goto L205;
    }
/*  fin ajoute FH jan 2004  */


  i_1 = *nbs + nbsd + nbsd;
  for (i = 1; i <= i_1; ++i)
     {
       w[i] = 0;
     }
  i_1 = nbsd;
  for (i = 1; i <= i_1; ++i)
     {
       a = (i_2 = sd[(i << 1) + 1], aabs (i_2));
/* Computing MIN */
       i_2 = arete[(a << 1) + 1], i_3 = arete[(a << 1) + 2];
       s1 = amin (i_2, i_3);
       w[i + i] = w[s1 + nbsd + nbsd];
       w[s1 + nbsd + nbsd] = i;
     }
  i_1 = nbt;
  for (t = 1; t <= i_1; ++t)
     {
       reft[t] = nothing;
       if (nu[t * 6 + 6] != 0)
	  {
	    for (i = 1; i <= 3; ++i)
	       {
/* Computing MIN */
		 i_2 = nu[i + t * 6], i_3 = nu[p3[i - 1] + t * 6];
		 ss1 = amin (i_2, i_3);
		 jsd = nbsd + nbsd + ss1;
	       L160:
		 isd = w[jsd];
		 if (isd > 0)
		    {
		      a = sd[(isd << 1) + 1];
		      if (a > 0)
			 {
			   if (nu[i + t * 6] == arete[(a << 1) + 1] && nu[p3[i -
					 1] + t * 6] == arete[(a << 1) + 2])
			      {
				reft[t] = sd[(isd << 1) + 2];
				w[isd + isd - 1] = t;
				w[jsd] = w[isd + isd];
				goto L170;
			      }
			 }
		      else if (a < 0)
			 {
			   if (nu[i + t * 6] == arete[(-a << 1) + 2] && nu[p3[i
				      - 1] + t * 6] == arete[(-a << 1) + 1])
			      {
				reft[t] = sd[(isd << 1) + 2];
				w[isd + isd - 1] = t;
				w[jsd] = w[isd + isd];
				goto L170;
			      }
			 }
		      else
			 {
			   *err = 11;
			 }
		      jsd = isd + isd;
		      goto L160;
		    }
	       L170:
		 ;
	       }
	  }
     }
  i_1 = nbsd;
  for (isd = 1; isd <= i_1; ++isd)
     {
       if (w[isd + isd - 1] == 0)
	  {
	    *err = 11;
	  }
       else
	  {
	    w[isd + isd] = 3;
	  }
     }
  if (*err != 0)
     {
       return 0;
     }
  i = nbsd + nbsd;
L200:
  if (i > 0)
     {
       ++w[i];
       if (w[i] <= 6)
	  {
	    ta = nu[w[i] + w[i - 1] * 6];
	    if (ta > 0)
	       {
		 ta /= 8;
		 if (nu[ta * 6 + 1] > 0)
		    {
		      nu[ta * 6 + 1] = -nu[ta * 6 + 1];
		      if (reft[ta] != reft[w[i - 1]])
			 {
			   if (reft[ta] != nothing)
			      {
			      }
			   else
			      {
				reft[ta] = reft[w[i - 1]];
			      }
			   w[i + 1] = ta;
			   w[i + 2] = 3;
			   i += 2;
			 }
		    }
	       }
	  }
       else
	  {
	    i += -2;
	  }
       goto L200;
     }
L205:     
  i_1 = nbt;
  for (ie = 1; ie <= i_1; ++ie)
     {
       if (nu[ie * 6 + 1] < 0)
	  {
	    nu[ie * 6 + 1] = -nu[ie * 6 + 1];
	  }
       else
	  {
	    for (i = 1; i <= 6; ++i)
	       {
		 nu[i + ie * 6] = 0;
	       }
	  }
     }
  return 1;
}				/* mshfrt_ */

/* ********************************************************************** */
int 
mshfr1_ (long *c, long *nu, long *it1, long *ita, long *is1, long *s2, long *err)
{
  /* Initialized data */

  static long     p3[5] =
  {2, 3, 1, 2, 3};

  static long     nbac, t, x, y, l1, l2, l3, s1, s3;

//    extern int mshfr2_();
  static long     la, ta;
  static long     s2t, s3t, det, lst[768] /* was [3][256] */ ;

  /* Parameter adjustments */
  nu -= 7;
  c -= 3;

  /* Function Body */
  t = *it1;
  s1 = nu[*is1 + t * 6];
  x = c[(*s2 << 1) + 1] - c[(s1 << 1) + 1];
  y = c[(*s2 << 1) + 2] - c[(s1 << 1) + 2];
  nbac = 0;
  l1 = *is1;
  l2 = p3[l1 - 1];
  l3 = p3[l2 - 1];
  s2t = nu[l2 + t * 6];
  s3t = nu[l3 + t * 6];
  la = l2 + 3;
L20:
  ++nbac;
  if (nbac > 256)
     {
       *err = 8;
       return 0;
     }
  lst[nbac * 3 - 2] = t;
  lst[nbac * 3 - 1] = la;
  ta = nu[la + t * 6];
  if (ta <= 0)
     {
       *err = 9;
       return 0;
     }
  t = ta / 8;
  la = ta - (t << 3);
  s3 = nu[p3[la - 3] + t * 6];
  if (s3 != *s2)
     {
       det = x * (c[(s3 << 1) + 2] - c[(s1 << 1) + 2]) - y * (c[(s3 << 1) +
						     1] - c[(s1 << 1) + 1]);
       if (det > 0)
	  {
	    la = p3[la - 4] + 3;
	  }
       else if (det < 0)
	  {
	    la = p3[la - 3] + 3;
	  }
       else
	  {
	    *err = 10;
	    return 0;
	  }
       goto L20;
     }
  mshfr2_ (&c[3], &nu[7], lst, &nbac, it1, ita, &s1, s2, err);
  return 0;
}				/* mshfr1_ */

/* ********************************************************************** */
int 
mshfr2_ (long *c, long *nu, long *lst, long *nbac, long *t, long *ta,
	 long *ss1, long *ss2, long *err)
{
  /* Initialized data */

  static long     mod3[3] =
  {2, 3, 1};

  /* System generated locals */
  long            i_1;

  /* Local variables */
  static long     i, x, y, a1, a2, pplst, s1, pslst, ptlst, s2, s3, s4,
                  ttlst, t1, t2, aa, i11, i12, i13, i21, i22, i23, x41,
                  y41, tt;

//    extern int mshopt_();
  static long     tt1, aas, det1, det2, det3, det4;

  /* Parameter adjustments */
  lst -= 4;
  nu -= 7;
  c -= 3;

  /* Function Body */
  x = c[(*ss1 << 1) + 1] - c[(*ss2 << 1) + 1];
  y = c[(*ss1 << 1) + 2] - c[(*ss2 << 1) + 2];
  i_1 = *nbac - 1;
  for (i = 1; i <= i_1; ++i)
     {
       lst[i * 3 + 1] = i + 1;
     }
  lst[*nbac * 3 + 1] = 0;
  ttlst = 1;
L20:
  ptlst = ttlst;
  pplst = 0;
L30:
  if (ptlst > 0)
     {
       t1 = lst[ptlst * 3 + 2];
       a1 = lst[ptlst * 3 + 3];
       tt1 = nu[a1 + t1 * 6];
       t2 = tt1 / 8;
       a2 = tt1 - (t2 << 3);
       i11 = a1 - 3;
       i12 = mod3[i11 - 1];
       i13 = mod3[i12 - 1];
       i21 = a2 - 3;
       i22 = mod3[i21 - 1];
       i23 = mod3[i22 - 1];
       s1 = nu[i13 + t1 * 6];
       s2 = nu[i11 + t1 * 6];
       s3 = nu[i12 + t1 * 6];
       s4 = nu[i23 + t2 * 6];
       x41 = c[(s4 << 1) + 1] - c[(s1 << 1) + 1];
       y41 = c[(s4 << 1) + 2] - c[(s1 << 1) + 2];
       det2 = (c[(s2 << 1) + 1] - c[(s1 << 1) + 1]) * y41 - (c[(s2 << 1) + 2]
						  - c[(s1 << 1) + 2]) * x41;
       det3 = (c[(s3 << 1) + 1] - c[(s1 << 1) + 1]) * y41 - (c[(s3 << 1) + 2]
						  - c[(s1 << 1) + 2]) * x41;
       if (det2 > 0 && det3 < 0)
	  {
/*         le quadrilataire est convexe on le retourne */
/*         update des sommets */
/* ------------------------- */
	    nu[i12 + t1 * 6] = s4;
	    nu[i22 + t2 * 6] = s1;
/*         update du pointeur suivant */
/* ----------------------------------- */
	    pslst = lst[ptlst * 3 + 1];
	    if (pslst > 0)
	       {
		 aas = lst[pslst * 3 + 3];
		 if (aas == i22 + 3)
		    {
		      lst[pslst * 3 + 2] = t1;
		      lst[pslst * 3 + 3] = i11 + 3;
		    }
	       }
/*         update des aretes a1,a2 */
/* ------------------------------- */
	    tt1 = nu[i22 + 3 + t2 * 6];
	    nu[a1 + t1 * 6] = tt1;
	    if (tt1 > 0)
	       {
		 tt = tt1 / 8;
		 aa = tt1 - (tt << 3);
		 nu[aa + tt * 6] = a1 + (t1 << 3);
	       }
	    else if (tt1 != nothing)
	       {
		 nu[-tt1 * 6 + 2] = a1 + (t1 << 3);
	       }
	    tt1 = nu[i12 + 3 + t1 * 6];
	    nu[a2 + t2 * 6] = tt1;
	    if (tt1 > 0)
	       {
		 tt = tt1 / 8;
		 aa = tt1 - (tt << 3);
		 nu[aa + tt * 6] = a2 + (t2 << 3);
	       }
	    else if (tt1 != nothing)
	       {
		 nu[-tt1 * 6 + 2] = a2 + (t2 << 3);
	       }
	    nu[i12 + 3 + t1 * 6] = i22 + 3 + (t2 << 3);
	    nu[i22 + 3 + t2 * 6] = i12 + 3 + (t1 << 3);
	    det1 = (c[(s1 << 1) + 1] - c[(*ss1 << 1) + 1]) * y - (c[(s1 << 1)
					     + 2] - c[(*ss1 << 1) + 2]) * x;
	    det4 = (c[(s4 << 1) + 1] - c[(*ss1 << 1) + 1]) * y - (c[(s4 << 1)
					     + 2] - c[(*ss1 << 1) + 2]) * x;
	    if (det1 < 0 && det4 > 0)
	       {
/*           le sommets s4 est dans omega */
		 lst[ptlst * 3 + 2] = t2;
		 lst[ptlst * 3 + 3] = i22 + 3;
	       }
	    else if (det1 > 0 && det4 < 0)
	       {
/*           le sommets s1 est dans omega */
		 lst[ptlst * 3 + 2] = t1;
		 lst[ptlst * 3 + 3] = i12 + 3;
	       }
	    else
	       {
		 if (pplst == 0)
		    {
		      ttlst = lst[ptlst * 3 + 1];
		      ptlst = ttlst;
		    }
		 else
		    {
		      ptlst = lst[ptlst * 3 + 1];
		      lst[pplst * 3 + 1] = ptlst;
		    }
		 goto L30;
	       }
	  }
       pplst = ptlst;
       ptlst = lst[ptlst * 3 + 1];
       goto L30;
     }
  if (ttlst != 0)
     {
       goto L20;
     }
  nu[i12 + 3 + t1 * 6] = nothing;
  nu[i22 + 3 + t2 * 6] = nothing;
  *t = t2;
  *ta = t1;
  i_1 = *nbac;
  for (i = 1; i <= i_1; ++i)
     {
       mshopt_ (&c[3], &nu[7], &lst[i * 3 + 2], 4L, err);
       mshopt_ (&c[3], &nu[7], &lst[i * 3 + 2], 5L, err);
       mshopt_ (&c[3], &nu[7], &lst[i * 3 + 2], 6L, err);
     }
  return 1;
}				/* mshfr2_ */


}
