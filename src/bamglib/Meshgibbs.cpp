//File gibbs2.cpp
/**************DO NOR REMOVE THIS BANNER***************/
/*  FreeFEM : Language for a Finite Element Method    */
/*  -------    Release 1.0:  June 1994.               */
/*  Authors: D. Bernardi, Y. Darmaillac F. Hecht,     */
/*           O. Pironneau                             */
/*  You may copy freely these files and use it for    */
/* teaching or research. These or part of these may   */
/* not be sold or used for a commercial purpose with- */
/* out our consent : fax (33)1 44 27 44 11            */
/* (e-mail)    Olivier.Pironneau@ann.jussieu.fr       */
/******************************************************/


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "vect.h"

#define mmax(a,b)(a>b?a:b)
#define mmin(a,b)(a<b?a:b)
#define ffalse 0
#define ttrue 1

/*  -- translated by f2c (version of 23 May 1992  14:18:33).
   You must link the resulting object file with the libraries:
	-lF77 -lI77 -lm -lc   (in that order)
*/


#define integer long
#define logical long

int gibbs1_(integer* n,integer*  record,integer*  ptvois);
int gibbs2_(integer* n,integer*  record,integer*  criter);
int gibbsa_(integer* n,integer*  ptvois,integer*  vois,integer*  r,integer*  m,
             integer*  nv,integer*  nx,integer*  ny,integer*  nn,integer*  w1,integer*  w2, 
	integer* pfold,integer*  pfnew,integer*  impre,integer*  nfout);
int gibbsb_(integer* x,integer*  y,integer*  n,integer*  ptvois,
						integer*  vois,integer*  nx,integer*  ny,integer*  nv,integer*  nn,integer*  m,
						integer*  wh,integer*  wl,integer* r, integer* impre, integer* nfout);
int gibbsc_(integer* nz,integer*  nv,integer*  niveau,integer*  n,integer* );
int gibbsd_(integer* racine,integer*  n,integer*  ptvois,integer*  
							vois,integer*  nv,integer*  r,integer*  niveau);
int gibbst_(integer* n,integer*  p,integer*  nv,integer*  nn,integer*  ptvois,integer*  vois,
						integer*  m,integer*  r,integer*  new_,integer*  option, 
						integer* pfnew,integer*  impre,integer*  nfout);

/* Subroutine */ int gibbs1_(integer* n,integer*  record,integer*  ptvois)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer crit, i, j, l, r, rec;

/* -----------------------------------------------------------------------
 */
/*     routine appele par gibbs0 */
/* -----------------------------------------------------------------------
 */
/*     but: trie record (ensemble de n sommet) telle que l'ordre des somm 
*/
/*     soit croissant (ordre du sommet i est ptvois(i+1)-ptvois(i)) */
/* -----------------------------------------------------------------------
 */

    /* Parameter adjustments */
    --ptvois;
    --record;

    /* Function Body */
    if (*n <= 1) {
	return 0;
    }
    l = *n / 2 + 1;
    r = *n;
L2:
    if (l <= 1) {
	goto L20;
    }
    --l;
    rec = record[l];
    crit = ptvois[record[l] + 1] - ptvois[record[l]];
    goto L3;
L20:
    rec = record[r];
    crit = ptvois[record[r] + 1] - ptvois[record[r]];
    record[r] = record[1];
    --r;
    if (r == 1) {
	goto L999;
    }
L3:
    j = l;
L4:
    i = j;
    j <<= 1;
    if ((i__1 = j - r) < 0) {
	goto L5;
    } else if (i__1 == 0) {
	goto L6;
    } else {
	goto L8;
    }
L5:
    if (ptvois[record[j] + 1] - ptvois[record[j]] < ptvois[record[j + 1] + 1] 
	    - ptvois[record[j + 1]]) {
	++j;
    }
L6:
    if (crit >= ptvois[record[j] + 1] - ptvois[record[j]]) {
	goto L8;
    }
    record[i] = record[j];
    goto L4;
L8:
    record[i] = rec;
    goto L2;
L999:
    record[1] = rec;
    return 0;
} /* gibbs1_ */

/* Subroutine */ int gibbs2_(integer* n,integer*  record,integer*  criter)
{
    static integer crit, i, j, l, r, rec;


/*     trie record selon les valeurs de criter(record(.)) croissantes */


    /* Parameter adjustments */
    --criter;
    --record;

    /* Function Body */
    if (*n <= 1) {
	return 0;
    }
    l = *n / 2 + 1;
    r = *n;
L2:
    if (l <= 1) {
	goto L20;
    }
    --l;
    rec = record[l];
    crit = criter[rec];
    goto L3;
L20:
    rec = record[r];
    crit = criter[rec];
    record[r] = record[1];
    --r;
    if (r == 1) {
	goto L999;
    }
L3:
    j = l;
L4:
    i = j;
    j <<= 1;
    if (j - r < 0) {
	goto L5;
    } else if (j == r) {
	goto L6;
    } else {
	goto L8;
    }
L5:
    if (criter[record[j]] < criter[record[j + 1]]) {
	++j;
    }
L6:
    if (crit >= criter[record[j]]) {
	goto L8;
    }
    record[i] = record[j];
    goto L4;
L8:
    record[i] = rec;
    goto L2;
L999:
    record[1] = rec;
    return 0;
} /* gibbs2_ */

/* Subroutine */ 
int gibbsa_(integer* n,integer*  ptvois,integer*  vois,integer*  r,integer*  m,
             integer*  nv,integer*  nx,integer*  ny,integer*  nn,integer*  w1,integer*  w2, 
	integer* pfold,integer*  pfnew,integer*  impre,integer*  nfout)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;

    /* Builtin functions */
    /* Subroutine */ int s_stop();

    /* Local variables */
    static integer nbcc, degi, bold, bnew, i, j, k, p, degre, x, y, p1, p2;
/*    extern  Subroutine  int gibbs1_();*/
    static integer pf;
/*    extern  Subroutine int gibbsb_(), gibbsd_(), gibbst_();*/
    static integer nbpass, niveau, pf1, option, old, new_, opt, new1;

/* -----------------------------------------------------------------------
 */
/*  but: calculer une renumerotation des sommets d'un graphe defini par: 
*/
/*     par la methode de gibbs */
/* -----------------------------------------------------------------------
 */
/*  entree */
/* -------- */
/*     n = nb de sommet du graphe */
/*      les voisins d'un sommet i ont pour numero : */
/*     ( vois(j) , j=ptvois(i),ptvois(i+1)-1 ) */

/*     impre   parametre d'impression */
/*     nfout   numero du fichier pour impression */

/*  sortie */
/*  ------ */
/*     r(1:n) tableau donnant la nouvelle numerotation: */
/*       r(i) = nouveau numero du sommet i */
/*     pfolf = ancien  profile */
/*     pfnew = nouveau profile */

/*  tableau de travail : */
/*  -------------------- */
/*     m(n) */
/*     nv(0:n+n) */
/*     nx(n) */
/*     ny(n) */
/*     nn(0:n) */
/*     w1(n) */
/*     w2(n) */

/* -----------------------------------------------------------------------
 */
/*     programmeur f. hecht  le 3/02/1987 */
/* -----------------------------------------------------------------------
 */

/*     tri des voisins d'un sommet du graphe par degre croissant */
/* --------------------------------------------------------------- */
    /* Parameter adjustments */
    --w2;
    --w1;
    --ny;
    --nx;
    --m;
    --r;
    --vois;
    --ptvois;

    /* Function Body */
    p2 = ptvois[1] - 1;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	p1 = p2 + 1;
	p2 = ptvois[i + 1] - 1;
	i__2 = p2 - p1 + 1;
	gibbs1_(&i__2, &vois[p1], &ptvois[1]);
/*       if(impre.le.-9) then */
/*        write (nfout,*) 'les voisin de ',i,'sont: ', (vois(j),j=p1,p
2) */
/*       endif */
/* L10: */
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	r[i] = 0;
/* L20: */
    }
/*     boucle sur les composante connexe du graphe */
    new_ = 0;
    nbcc = 0;
L30:
    if (new_ < *n) {
	++nbcc;
/*       recherche d'une racine y (un sommet non numerote) de degree m
ini */
	y = 0;
	degre = *n + 1;
	i__1 = *n;
	for (i = 1; i <= i__1; ++i) {
	    if (r[i] <= 0) {
		degi = ptvois[i + 1] - ptvois[i];
		if (degi < degre) {
		    degre = degi;
		    y = i;
		}
	    }
/* L40: */
	}
	if (y == 0) {
	   return -3;/*  s_stop("fatal erreur  gibbs 2 : pb racine", 33L); */
	}
	gibbsd_(&y, n, &ptvois[1], &vois[1], nv, &r[1], &niveau);
	nbpass = 0;
L50:
	++nbpass;
	x = y;
	p = niveau;
	k = 0;
	i__1 = nv[p + 1];
	for (i = nv[p] + 1; i <= i__1; ++i) {
	    ++k;
	    m[k] = nv[i];
/* L60: */
	}
	gibbs1_(&k, &m[1], &ptvois[1]);
	i__1 = k;
	for (i = 1; i <= i__1; ++i) {
	    y = m[i];
	    gibbsd_(&y, n, &ptvois[1], &vois[1], nv, &r[1], &niveau);
	    if (niveau > p) {
		goto L50;
	    }
/* L70: */
	}
	y = m[1];
/*        if(impre.lt.0) then */
/*          write(nfout,*) */
/*     +  '    nb de pass pour trouver le pseudo diametre',nbpass */
/*     +         ,' x=',x,',y=',y,' de la composante connexe ',nbcc */

/*          write (nfout,*) ('-',i=1,78) */
/*        endif */
/*       optimisation de la descendance de la numerotation */
/*       ------------------------------------------------- */
	gibbsb_(&x, &y, n, &ptvois[1], &vois[1], &nx[1], &ny[1], nv, nn, &m[1]
		, &w1[1], &w2[1], &r[1], impre, nfout);

/*     renumerotation de cuthill mac kee avec la meilleur des 4 option
s */
/*     --------------------------------------------------------------
--- */
	pf = 1073741824;
	option = -2;
	new1 = new_;
	for (opt = -2; opt <= 2; ++opt) {
	    new_ = new1;
	    if (opt != 0) {
		gibbst_(n, &p, nv, nn, &ptvois[1], &vois[1], &m[1], &r[1], &
			new_, &opt, &pf1, impre, nfout);
		if (pf1 < pf) {
		    pf = pf1;
		    option = opt;
		}
	    }
/* L80: */
	}
/*        if(impre.ne.0) write (nfout,*) '    on a choisi l''option ',
 */
/*     +             option,', new =',new */
	new_ = new1;
	gibbst_(n, &p, nv, nn, &ptvois[1], &vois[1], &m[1], &r[1], &new_, &
		option, &pf1, impre, nfout);
	goto L30;
    }
/*      if(impre.ne.0) write(nfout,*) */
/*     +       '   nb de composante connexe du graphe =',nbcc */
/*     calcul du profile */
    *pfold = 0;
    *pfnew = 0;
    bnew = 0;
    bold = 0;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	old = i;
	new_ = r[i];
	i__2 = ptvois[i + 1] - 1;
	for (j = ptvois[i]; j <= i__2; ++j) {
/* Computing MIN */
	    i__3 = old, i__4 = vois[j];
	    old = mmin(i__3,i__4);
/* Computing MIN */
	    i__3 = new_, i__4 = r[vois[j]];
	    new_ = mmin(i__3,i__4);
/* L100: */
	}
	*pfold = *pfold + i - old + 1;
/* Computing MAX */
	i__2 = bold, i__3 = i - old + 1;
	bold = mmax(i__2,i__3);
	*pfnew = *pfnew + r[i] - new_ + 1;
/* Computing MAX */
	i__2 = bnew, i__3 = r[i] - new_ + 1;
	bnew = mmax(i__2,i__3);
/* L110: */
    }
/*      if(impre.ne.0) then */
/*        write(nfout,*)'profile  old  = ',pfold,', profile  new = ',pfnew
 */
/*        write(nfout,*)'1/2 bande old = ',bold ,', 1/2 band new = ',bnew 
*/
/*      endif */
return 0;
} /* gibbsa_ */

/* Subroutine */ int gibbsb_(integer* x,integer*  y,integer*  n,integer*  ptvois,
integer*  vois,integer*  nx,integer*  ny,integer*  nv,integer*  nn,integer*  m,
integer*  wh,integer*  wl,integer* r, integer* impre, integer* nfout)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static  flag_;
    static integer i, j, k, p, s, h0, i1, l0, i2;
/*    extern  Subroutine  int gibbs1_(); */
    static integer lg;
/*    extern  Subroutine  int gibbsd_(), gibbsc_();*/
    static integer niveau, mxcanx, mxcany, nbc;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 */
/* ...................................................................... 
*/
/*     attention on met la descente optimiser dans r <0 ou nulle */
/* .......................................................................
 */
    /* Parameter adjustments */
    --r;
    --m;
    --ny;
    --nx;
    --vois;
    --ptvois;

    /* Function Body */
    gibbsd_(y, n, &ptvois[1], &vois[1], nv, &r[1], &niveau);
    gibbsc_(&ny[1], nv, &niveau, n, &mxcany);
    gibbsd_(x, n, &ptvois[1], &vois[1], nv, &r[1], &niveau);
    p = niveau;
    gibbsc_(&nx[1], nv, &niveau, n, &mxcanx);
    flag_ = ffalse;
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	if (nx[i] + ny[i] == p) {
	    r[i] = -nx[i];
	} else if (nx[i] >= 0) {
	    flag_ = ttrue;
	    r[i] = -1073741824;
	} else {
	    if (r[i] <= 0) {
		r[i] = -1073741822;
	    }
	}
/* L20: */
    }
    if (flag_) {
/*       calcul des composantes connexe du graphe sans les sommets de 
nn */
/*       ------------------------------------------------------------
--- */
	j = *n;
	k = 0;
	nbc = 0;
	nv[nbc] = j;
L30:
	++k;
	if (k <= *n) {
	    if (r[k] == -1073741824) {
/*           recherche de la fermeture transitive partant de k
 */
		++nbc;
		i = -1;
		s = k;
L40:
		++i;
		wl[i] = ptvois[s];
		wh[i] = ptvois[s + 1];
		++j;
		nv[j] = s;
		r[s] = -1073741823;
L50:
		if (i >= 0) {
		    if (wl[i] < wh[i]) {
			s = vois[wl[i]];
			++wl[i];
			if (r[s] == -1073741824) {
			    goto L40;
			}
			goto L50;
		    }
		    --i;
		    goto L50;
		}
		nv[nbc] = j;
		m[nbc] = nbc;
	    }
	    goto L30;
	}
/*        if(impre.lt.0) write(nfout,*) */
/*     +         ' nb de composante connexe du graphe reduit =',nbc */

/* --------------- fin de construction des composantes connexes------
--- */
/*        nv(0)=n */
/*        if(impre.le.-10) write(nfout,5555)'nv(0:n+n) = ',(nv(i),i=0,
n+n) */
	gibbs1_(&nbc, &m[1], nv);
/*        if(impre.le.-10)write(nfout,5555)'trie m =',(m(i),i=1,nbc) 
*/
	i__1 = p;
	for (i = 0; i <= i__1; ++i) {
	    nn[i] = 0;
/* L60: */
	}
	i__1 = *n;
	for (i = 1; i <= i__1; ++i) {
	    j = -r[i];
	    if (j >= 0 && j <= p) {
		++nn[j];
	    }
/* L70: */
	}

/*       boucle sur les composante connexes par ordre croissantes */
/*       -------------------------------------------------------- */
	for (k = nbc; k >= 1; --k) {
	    i = m[k];
	    i1 = nv[i - 1] + 1;
	    i2 = nv[i];
	    lg = i2 - i1 + 1;
/*         if(impre.le.-7) */
/*     +       write(nfout,*) k,' composante ',i,',lg=',lg,',i1,i2
=',i1,i2 */
/*         if(impre.le.-8) */
/*     +       write (nfout,5555)' ',(nv(i),i=i1,i2) */
	    h0 = 0;
	    l0 = 0;
	    i__1 = p;
	    for (j = 0; j <= i__1; ++j) {
		wh[j] = nn[j];
		wl[j] = nn[j];
/* L90: */
	    }
	    i__1 = i2;
	    for (i = i1; i <= i__1; ++i) {
		s = nv[i];
		++wh[nx[s]];
		++wl[p - ny[s]];
/* L100: */
	    }
	    i__1 = p;
	    for (j = 0; j <= i__1; ++j) {
		if (wh[j] != nn[j]) {
/* Computing MAX */
		    i__2 = wh[j];
		    h0 = mmax(i__2,h0);
		}
		if (wl[j] != nn[j]) {
/* Computing MAX */
		    i__2 = wl[j];
		    l0 = mmax(i__2,l0);
		}
/* L110: */
	    }
	    if (h0 < l0 || h0 == l0 && mxcanx <= mxcany) {
/*           if(impre.le.-2) write(nfout,*) */
/*     +       '         h0 = ',h0,',l0 = ',l0,'  ------- XXXX
 --------' */
		i__1 = i2;
		for (i = i1; i <= i__1; ++i) {
		    s = nv[i];
		    r[s] = -nx[s];
		    ++nn[-r[s]];
/* L120: */
		}
	    } else {
/*           if (impre.le.-2) write(nfout,*) */
/*     +       '         h0 = ',h0,',l0 = ',l0,'  ------- YYYY
 --------' */
		i__1 = i2;
		for (i = i1; i <= i__1; ++i) {
		    s = nv[i];
		    r[s] = -p + ny[s];
		    ++nn[-r[s]];
/* L130: */
		}
	    }
/* L140: */
	}
    }
/*     on met les nouveaux niveaux de la descendance optimiser dans nn */
/*     ----------------------------------------------------------------- 
*/
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	if (r[i] > 0) {
	    nn[i] = -1;
	} else if (r[i] == -1073741822) {
	    nn[i] = -2;
	} else {
	    nn[i] = -r[i];
	}
/* L150: */
    }
/*      if(impre.le.-10)write (nfout,5555)' nn(i)=',(nn(i),i=1,n) */
/* 5555  format('            --------   ',a,/,5(15x,10(i5)/)) */
return 0;} /* gibbsb_ */

/* Subroutine */ int gibbsc_(integer* nz,integer*  nv,integer*  niveau,integer*  n,integer*  mxz)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    static integer i, j;

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 */
    /* Parameter adjustments */
    --nz;

    /* Function Body */
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	nz[i] = -1;
/* L10: */
    }
    *mxz = 0;
    i__1 = *niveau;
    for (i = 0; i <= i__1; ++i) {
/* Computing MAX */
	i__2 = *mxz, i__3 = nv[i + 1] - nv[i];
	*mxz = mmax(i__2,i__3);
	i__2 = nv[i + 1];
	for (j = nv[i] + 1; j <= i__2; ++j) {
	    nz[nv[j]] = i;
/* L20: */
	}
    }
return 0;} /* gibbsc_ */

/* Subroutine */ int gibbsd_(integer* racine,integer*  n,integer*  ptvois,integer*  vois,integer*  nv,integer*  r,integer*  niveau)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i, k, s, sv, stk, stk1, stk2;

/* -----------------------------------------------------------------------
 */
/*     but construire la structure des descendant de racine  du graphe */
/* -----------------------------------------------------------------------
 */
/*     sortie : */
/*     -------- */
/*     nv est la structure des niveaux */
/*     les sommets du niveau (i =0,niveau_ sont defini par : */
/*        (nv(j),j=nv(i),nv(i+1)-1) */

/*     le tableau r(i) n'est modifier que sur les sommets */
/*       de la composante connexe du graphe contenant la racine */

/* -----------------------------------------------------------------------
 */

/*     on demark tout les sommets non remuneroter */
/* -------------------------------------------------- */
    /* Parameter adjustments */
    --r;
    --vois;
    --ptvois;

    /* Function Body */
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	if (r[i] < 0) {
	    r[i] = 0;
	}
/* L10: */
    }

/*    initialisation */

    stk = *n - 1;
    nv[0] = stk;
    stk2 = stk;
    *niveau = 0;
    ++stk;
    nv[stk] = *racine;
    r[*racine] = -1;
L20:
    if (stk2 < stk) {
	++(*niveau);
	stk1 = stk2 + 1;
	nv[*niveau] = stk;
	stk2 = stk;
/*        print *,' ------- niveau =',niveau,' stk=',stk1,stk2 */
	i__1 = stk2;
	for (k = stk1; k <= i__1; ++k) {
	    s = nv[k];
/*         print *,'----------------- s=',s */
	    i__2 = ptvois[s + 1] - 1;
	    for (i = ptvois[s]; i <= i__2; ++i) {
/*               pour tout les sommets (sv) voisin */
/*                d'un sommet (s) du niveau precedent */
		sv = vois[i];
/*          print *,' voisin =',sv */
/*               si le sommet n'est pas marque on le marque et
 on l'ajout */
		if (r[sv] == 0) {
		    ++stk;
		    nv[stk] = sv;
		    r[sv] = -1;
		}
/* L30: */
	    }
/* L40: */
	}
	goto L20;
    }
    --(*niveau);
/*      call pnv(' gibbsd ',n,nv,niveau) */
return 0;} /* gibbsd_ */


/* Subroutine */ int gibbst_(integer* n,integer*  p,integer*  nv,integer*  nn,integer*  ptvois,integer*  vois,integer*  m,integer*  r,integer*  new_,integer*  option, 
	integer* pfnew,integer*  impre,integer*  nfout)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;

    /* Local variables */
    static integer nbsc, bnew, knew, step, plus, i, j, k, s, debut, i1, i2;
/*    extern Subroutine int gibbs2_();*/
    static integer fin;


/*     construction de la stucture de niveau dans nv a partir de nn */
/*     ------------------------------------------------------------ */
    /* Parameter adjustments */
    --r;
    --m;
    --vois;
    --ptvois;

    /* Function Body */
    nv[0] = *n;
    i__1 = *p + 1;
    for (i = 1; i <= i__1; ++i) {
	nv[i] = 0;
/* L150: */
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	if (nn[i] >= 0) {
	    ++nv[nn[i] + 1];
	}
/* L160: */
    }
    i__1 = *p;
    for (i = 0; i <= i__1; ++i) {
	nv[i + 1] += nv[i];
/* L170: */
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	if (nn[i] >= 0) {
	    j = nn[i];
	    ++nv[j];
	    nv[nv[j]] = i;
	}
/* L180: */
    }
    for (i = *p; i >= 0; --i) {
	nv[i + 1] = nv[i];
/* L190: */
    }
    nv[0] = *n;
    nbsc = nv[*p + 1] - nv[0];
/*     --- fin de la construction ------------------------------------ */
    if (*option == -2) {
	i__1 = *impre - 1;
    }
    i__1 = *n;
    for (i = 1; i <= i__1; ++i) {
	m[i] = *n * 3 + ptvois[i + 1] - ptvois[i];
/* L10: */
    }
    if ((((int)*option) == 1)||(((int)*option) == -1)) {
	debut = 0;
	fin = *p;
	step = 1;
    } else {
	debut = *p;
	fin = 0;
	step = -1;
    }
    i__1 = fin;
    i__2 = step;
    for (i = debut; i__2 < 0 ? i >= i__1 : i <= i__1; i += i__2) {
	i1 = nv[i] + 1;
	i2 = nv[i + 1];
	i__3 = i2 - i1 + 1;
	gibbs2_(&i__3, &nv[i1], &m[1]);
	i__3 = i2;
	for (j = i1; j <= i__3; ++j) {
	    s = nv[j];
	    i__4 = ptvois[s + 1] - 1;
	    for (k = ptvois[s]; k <= i__4; ++k) {
/* Computing MIN */
		i__5 = m[vois[k]];
		m[vois[k]] = mmin(i__5,j);
/* L20: */
	    }
/* L30: */
	}
/* L40: */
    }
    if (*option > 0) {
	knew = *new_;
	plus = 1;
    } else {
	knew = *new_ + nbsc + 1;
	plus = -1;
    }
    *new_ += nbsc;
/*      if(option.gt.0) then */
/*        do 60 k = debut , fin , step */
/*          do 60 j = nv(k+1),nv(k)+1,-1 */
/*            knew = knew + plus */
/*            r(nv(j)) = knew */
/* 60      continue */
/*      else */
    i__2 = fin;
    i__1 = step;
    for (k = debut; i__1 < 0 ? k >= i__2 : k <= i__2; k += i__1) {
	i__3 = nv[k + 1];
	for (j = nv[k] + 1; j <= i__3; ++j) {
	    knew += plus;
	    r[nv[j]] = knew;
/* L70: */
	}
    }
/*      endif */
    *pfnew = 0;
    bnew = 0;
    i__3 = *n;
    for (i = 1; i <= i__3; ++i) {
	k = r[i];
	if (k > 0) {
	    i__1 = ptvois[i + 1] - 1;
	    for (j = ptvois[i]; j <= i__1; ++j) {
		if (r[vois[j]] > 0) {
/* Computing MIN */
		    i__2 = k, i__4 = r[vois[j]];
		    k = mmin(i__2,i__4);
		}
/* L100: */
	    }
	    *pfnew = *pfnew + r[i] - k + 1;
/* Computing MAX */
	    i__1 = bnew, i__2 = r[i] - k + 1;
	    bnew = mmax(i__1,i__2);
	}
/* L110: */
    }
/*      if(impre.lt.0.or.impre.gt.2) then */
/*        write(nfout,*) '      option =',option,', profile =',pfnew */
/*     +       ,', 1/2 bande =',bnew,', new=',new,', nbss composante=',nbsc
 */
/*      endif */
return 0;} /* gibbst_ */

/* function */ 
int Triangles::gibbsv (integer* ptvoi,
		       integer* vois,integer* lvois,integer* w,integer* v)
{
  /* System generated locals */
  integer  i__2;
  
  /* Local variables */
  integer i, j, k, T, ss, iii, ptv, ptv1;
  integer nbss = nbv; 
  /*--- Prepare les donees pour gibbsa en construisant ptvoi, vois, lvois -
    ------------*/
  /*     in */
  /*     ---   nbnt =3 pour des triangles 2D, */
  /* 			nbt =  nb de triangle */
  /*    		nbss = nb de sommets */
  /*           nsea = numeros de 3 sommets de chaque triangle (me) */
  /*     out */
  /*     --- 	ptvoi, vois, lvois, err */
  /*      tableaux de travail w, v */
  /*-----------------------------------------------------------------------
    ----------*/
  /* Parameter adjustments */
  --v;
  --w;
  --vois;
  --ptvoi;
  long nt = nbt-NbOutT;
  /* Function Body */
  for (i = 1; i <= nbss; ++i) {
    w[i] = -1;
    ptvoi[i] = 0; }
  ptvoi[nbss + 1] = 0;
  for (i = 0; i < nt; ++i) 
    { 
      assert(triangles[i].link);
      for (j = 0; j < 3; ++j)
	{
	  ss = Number(triangles[i][j])+1;
	  ++ptvoi[ss + 1];
	  w[ss] = 0;
	}
    }
  
  for (i = 1; i <= nbss; ++i) 
    ptvoi[i + 1] += ptvoi[i];
  
  for (i = 0; i < nt; ++i)
    if (triangles[i].link) 
      for (j = 0; j < 3; ++j) 
	{
	  ss = Number(triangles[i][j])+1;
	  ++ptvoi[ss];
	  v[ptvoi[ss]] = i;
	}
  
    ptv1 = 0;
    iii = 1;
    for (i = 1; i <= nbss; ++i) {
      ptv = ptv1 + 1;
      ptv1 = ptvoi[i];
      ptvoi[i] = iii;
      i__2 = ptv1;
      for (j = ptv; j <= i__2; ++j) {
	T = v[j];
	for (k = 0; k < 3; ++k) {
	  ss = Number(triangles[T][k])+1;  /*  nsea[k + T * nsea_dim1]; */
	  if (w[ss] != i) {
	    w[ss] = i;
	    if (iii > *lvois)  return 2 ;
	    /* print*,'pas assez de place memoire' */
	    
	    vois[iii] = ss;
	    ++iii;}
	}
      }
    }
    ptvoi[nbss + 1] = iii;
    *lvois = iii - 1;
    return 0; /* OK */
    return 0;} /* gibbsv_ */

int Triangles::gibbs()
/* -------- 
	renumber vertices by gibbs method; updates triangle and edge array
	in:   mesh  
	out:   mesh
 	auxiliary arrays: ptvois,vois,r,m,nv,nx,ny,nn,w1,w2,f 
 	all of size nv+1 except vois (10(nv+1)) and nv (2(nv+1))
 	err = -1 : memory alloc pb; err = -3: fatal erreur  gibbs 2 : pb racine
*/
{
  long nv = nbv;
  long nt = nbt-NbOutT;
    long i, j, pfold, pfnew;
    long* ptvois=NULL;
    long* vois=NULL;
    long* nn =NULL;
    long* r =NULL;
    long* m =NULL;
    long* nnv =NULL;
    long* nx =NULL;
    long* ny =NULL;
    long* w1 =NULL;
    long* w2=NULL;
    long nbvoisin =  10*nv;
    long printint=0, iodev=6;
    int err=0;
    ptvois = new long[nv+1]; 		//(long*)calloc((long)(nv + 1) , sizeof(long));
    nn = 	 new long[3*nt]; 			//(long*)calloc(3 * nt ,sizeof(long));
    vois = 	 new long[nbvoisin+10];	//(long*)calloc((long)(nbvoisin + 10) , sizeof(long)); 
    r = 	 new long[nv+1];				//(long*)calloc((long)(nv + 1) , sizeof(long));
    if((!ptvois)||(!nn)||(!vois)||(!r)) return -1;
    err = gibbsv(ptvois,vois,&nbvoisin,r,nn) ;
    delete [] nn;					// free(nn);
    if(err==0)
      {
	m = new long[nv+1];
	nn = new long[nv+1];
	nnv = new long[(nv+1)<<1];
	nx = new long[nv+1];
	ny = new long[nv+1];
	w1 = new long[nv+1];
	w2 = new long[nv+1];
	long lnv = nv;
	err = gibbsa_ (&lnv, ptvois, vois, r, m, nnv, nx, ny, nn, w1, w2, &pfold, &pfnew,
		      &printint, &iodev);
	delete [] m;
	delete [] nnv;
	delete [] nn;
	delete [] nx;
	delete [] ny;
	delete [] w1;
	delete [] w2;
      }
    
    delete [] vois;
  delete [] ptvois;
  /*          
	      if (err == 0 && (pfnew <= pfold))
	      {
	      A<bVertex> f(nv);
	      for (i = 0; i < nv; ++i)
	      {	f[i].x = v[i].x;
	      f[i].y = v[i].y;
	      f[i].where = v[i].where;
	      }
	      for (i = 0; i < nv; ++i)
	      {	v[r[i] - 1].x = f[i].x;
	 		v[r[i] - 1].y = f[i].y;
	 		v[r[i] - 1].where = f[i].where;
			}
			
       for (j = 0; j < nt; ++j)  // updates triangle array
       for (i = 0; i < 3; i++)
       t[j].v[i] = &v[r[no(t[j].v[i])] - 1];
       
       for (j = 0; j < ne; ++j)	// updates edge array
       {
	   		e[j].in = &v[r[no(e[j].in)] - 1];
	   		e[j].out = &v[r[no(e[j].out)] - 1];
			}
			f.destroy();
       if (!NumThinGrid) 
       {  NumThinGrid= new int [nv];
       for (i=0;i<nv;i++) NumThinGrid[i]=i;// Same numbering 
       }
       for (i=0;i<nv;i++) NumThinGrid[i]=r[NumThinGrid[i]]-1;  
       
       }
  */
  delete [] r;
  return err;
} 
			
/*  message d'erreur:         *err = 2;    print*,'pas assez de place memoire'   */
