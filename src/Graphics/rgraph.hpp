/**************DO NOR REMOVE THIS BANNER***************/

/*  FreeFEM : Language for a Finite Element Method    */

/*  -------    Release 1.0:  June 1994.               */

/*  Authors: D. Bernardi, Y. Darmaillac F. Hecht,     */

/*           O. Pironneau                             */

/*  You may copy freely these files and use it for    */

/* teaching or research. These or part of these may   */

/* not be sold or used for a commercial purpose with- */

/* out our consent : fax (33)1 44 27 44 11            */

/*  modified for bamg by F Hecht  dec 1997            */

/* add of  InPtScreen   and   InRecScreen             */

/* (e-mail)    Olivier.Pironneau@ann.jussieu.fr       */

/* (e-mail)    Frederic.Hecht@inria.fr                */

/******************************************************/

#ifndef RGRAPH_H_

#define RGRAPH_H_

#ifdef __cplusplus

#include "config.h"

extern "C" {

	void getcadre(float &xmin,float &xmax, float &ymin, float &ymax);

	void GetScreenSize(int &ix,int &iy);	// unused

	char Getxyc(float &x,float &y);

	void ShowHelp(const char * s,int k=1); //tlr: deplace ici

#endif

	void erreur(char *s);

	void initgraphique();	// only in main

	void closegraphique();

	void showgraphic();

	void rattente(int waitm);

	void cadre(float xmin,float xmax, float ymin, float ymax);

	void cadreortho(float centrex, float centrey, float rayon);

	void couleur(int c);

	int  LaCouleur();

	void pointe(float x, float y);

	int  InPtScreen(float x, float y);

	int  InRecScreen(float x1, float y1,float x2, float y2);

	void plotstring(const char *s);

	void rmoveto(float x, float y);

	void rlineto(float x, float y);

	void penthickness(int );

	int execute(const char* s);

	void reffecran();

	void fillpoly(int n, float *poly);

	void SetColorTable(int nb);

	float GetHeigthFont(); 

	void compile(char *fname);

	void compileString(char *texte);/*tlr: add a string stream */

	void openPS(const char * );

	void closePS(void);

	void coutmode(short i);

	void myexit(int err=0);

	void viderbuff();

	void Commentaire(const char *);

	void NoirEtBlanc(int NB);

	void MettreDansPostScript(int in);//  oui=1 ou non=0 

    int getprog(char* fn,int , char** argvptr);
    void setgrey(bool );
    int getgrey( );
    
	    

#ifdef __cplusplus

}



#ifdef TERM_USED

/** Ouput on the terminal window */	   

class myostream {

	int channel;

	public:

	myostream(int ch) { channel = ch; }

	myostream& operator<<(char c);

	myostream& operator<<(const char *s);

	myostream& operator<<(const void *p);

	myostream& operator<<(int n);

	myostream& operator<<(unsigned int n) { return *this << (int)n; }

	myostream& operator<<(long n);

	myostream& operator<<(unsigned long n) { return *this << (long)n; }

	myostream& operator<<(double n);

	myostream& operator<<(float n) { return *this << (double)n; }

	myostream& operator<<(__omanip func);

};



extern myostream termout; // could be cout, or another thing

extern myostream termerr; // could be cerr, or another thing

#define cout termout

#define cerr termerr

#endif /* TERM_USED */



#endif /* __cplusplus */



#endif /* RGRAPH_H_ */

