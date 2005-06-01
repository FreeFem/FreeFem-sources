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

// Modif F Hecht for dll on window 
// pass all the graphic function via pointeur
//  now the real graphic function of a pointeur xxxx is call  xxxx_
//   We just need to add INGRAPH in all graphic.cpp version
// ----- 
#ifdef INGRAPH
#define EXTERNFF(t,f,arg) t f##_ arg;t (*f) arg  = f##_; t f##_ arg 
#else
#define EXTERNFF(t,f,arg) extern t (*f) arg
#endif
#define EXTERN 
#ifdef __cplusplus

#include "config-wrapper.h"


//extern "C" {
EXTERNFF(	void ,getcadre,(float &xmin,float &xmax, float &ymin, float &ymax)) ;
EXTERNFF(	void ,GetScreenSize,(int &ix,int &iy)) ;	// unused
EXTERNFF(	char ,Getxyc,(float &x,float &y)) ;
EXTERNFF(	void ,ShowHelp,(const char * s,int k)) ; // k=1 ??//tlr: deplace ici
#endif
EXTERNFF(	void ,erreur,(char *s)) ;
EXTERNFF(	void ,initgraphique,()) ;	// only in main
EXTERNFF(	void ,closegraphique,()) ;
EXTERNFF(	void ,showgraphic,()) ;
EXTERNFF(	void ,rattente,(int waitm)) ;
EXTERNFF(	void ,cadre,(float xmin,float xmax, float ymin, float ymax)) ;
EXTERNFF(	void ,cadreortho,(float centrex, float centrey, float rayon)) ;
EXTERNFF(	void ,couleur,(int c)) ;
EXTERNFF(	int  ,LaCouleur,()) ;
EXTERNFF(	void ,pointe,(float x, float y)) ;
EXTERNFF(	int ,InPtScreen,(float x, float y)) ;
EXTERNFF(	int ,InRecScreen,(float x1, float y1,float x2, float y2)) ;
EXTERNFF(	void ,plotstring,(const char *s)) ;
EXTERNFF(	void ,rmoveto,(float x, float y)) ;
EXTERNFF(	void ,rlineto,(float x, float y)) ;
EXTERNFF(	void ,penthickness,(int )) ;
EXTERNFF(	int ,execute,(const char* s)) ;
EXTERNFF(	void ,reffecran,()) ;
EXTERNFF(	void ,fillpoly,(int n, float *poly)) ;
EXTERNFF(	void ,SetColorTable,(int nb)) ;
EXTERNFF(	float ,GetHeigthFont,()) ; 
//EXTERNFF(	void ,compile,(char *fname)) ;
//EXTERNFF(	void ,compileString,(char *texte)) ;/*tlr: add a string stream */
EXTERNFF(	void ,openPS,(const char * )) ;
EXTERNFF(	void ,closePS,(void)) ;
EXTERNFF(	void ,coutmode,(short i)) ;
EXTERNFF(	void ,myexit,(int err)) ; // err=0 ??
EXTERNFF(	void ,viderbuff,()) ;
EXTERNFF(	void ,Commentaire,(const char *)) ;
EXTERNFF(	void ,NoirEtBlanc,(int NB)) ;
EXTERNFF(	void ,MettreDansPostScript,(int in)) ;//  oui=1 ou non=0 
EXTERNFF(    int ,getprog,(char* fn,int , char** argvptr)) ;
EXTERNFF(    void ,setgrey,(bool )) ;
EXTERNFF(    int ,getgrey,(    )) ;


#ifdef  INGRAPH 
#define getcadre getcadre_
#define GetScreenSize GetScreenSize_
#define Getxyc Getxyc_
#define ShowHelp ShowHelp_
#define erreur erreur_
#define initgraphique initgraphique_
#define closegraphique closegraphique_
#define showgraphic showgraphic_
#define rattente rattente_
#define cadre cadre_
#define cadreortho cadreortho_
#define couleur couleur_
#define LaCouleur LaCouleur_
#define pointe pointe_
#define InPtScreen InPtScreen_
#define InRecScreen InRecScreen_
#define plotstring plotstring_
#define rmoveto rmoveto_
#define rlineto rlineto_
#define penthickness penthickness_
#define execute execute_
#define reffecran reffecran_
#define fillpoly fillpoly_
#define SetColorTable SetColorTable_
#define GetHeigthFont GetHeigthFont_
#define compile compile_
#define compileString compileString_
#define openPS openPS_
#define closePS closePS_
#define coutmode coutmode_
#define myexit myexit_
#define viderbuff viderbuff_
#define Commentaire Commentaire_
#define NoirEtBlanc NoirEtBlanc_
#define MettreDansPostScript MettreDansPostScript_
#define getprog getprog_
#define setgrey setgrey_
#define getgrey getgrey_
#endif
	    
#ifdef __cplusplus

//}



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

