// SUMMARY  :      
// USAGE    :        
// ORG      : 
// AUTHOR   : Frederic Hecht
// E-MAIL   : hecht@ann.jussieu.fr
/* -------------------                                */
/* (e-mail)    Olivier.Pironneau@ann.jussieu.fr       */
/* (e-mail)    hecht@ann.jussieu.fr                   */
/******************************************************/
/*
 
 This file is part of Freefem++
 
 Freefem++ is free software; you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation; either version 2.1 of the License, or
 (at your option) any later version.
 
 Freefem++  is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.
 
 You should have received a copy of the GNU Lesser General Public License
 along with Freefem++; if not, write to the Free Software
 Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */
// E-MAIL :   Frederic.Hecht@Inria.fr   
//
// ORIG-DATE:     Dec 97
#ifndef RGRAPH_H_

#define RGRAPH_H_

// Modif F Hecht for dll on win32 
// pass all the graphic function via pointeur
//  now the real graphic function of a pointeur xxxx is call  xxxx_
//   We just need to add FF_GRAPH_SET_PTR in all graphic.cpp version
//  the dcl of xxxx pointeur is  done via GRAPH_PTR_DCL macro (set one time ) 
//  
// ----- 
#ifdef FF_GRAPH_SET_PTR
#define EXTERNFF(t,f,arg) t f##_ arg;  extern t (*f) arg 
#else
#ifdef FF_GRAPH_PTR_DCL
#define EXTERNFF(t,f,arg)  t (*f) arg
#else
#define EXTERNFF(t,f,arg) extern t (*f) arg
#endif
#endif
//#define EXTERN 
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
EXTERNFF(	void ,SetColorTable1,(int nb,bool hsv,int nbcolors,float *colors)) ;
EXTERNFF(	float ,GetHeigthFont,()) ; 
//  old function for freefem+  
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


// wrapping of function  -----
#ifdef  FF_GRAPH_SET_PTR 
static int init_ff_graph_ptr_func()
{ //  a small function to set all pointeur
	getcadre=getcadre_;
	GetScreenSize=GetScreenSize_;
	Getxyc=Getxyc_;
	ShowHelp=ShowHelp_;
	erreur=erreur_;
	initgraphique=initgraphique_;
	closegraphique=closegraphique_;
	showgraphic=showgraphic_;
	rattente=rattente_;
	cadre=cadre_;
	cadreortho=cadreortho_;
	couleur=couleur_;
	LaCouleur=LaCouleur_;
	pointe=pointe_;
	InPtScreen=InPtScreen_;
	InRecScreen=InRecScreen_;
	plotstring=plotstring_;
	rmoveto=rmoveto_;
	rlineto=rlineto_;
	penthickness=penthickness_;
	execute=execute_;
	reffecran=reffecran_;
	fillpoly=fillpoly_;
	SetColorTable=SetColorTable_;
	SetColorTable1=SetColorTable1_;
	GetHeigthFont=GetHeigthFont_;
	//compile=compile_;
	//compileString=compileString_;
	openPS=openPS_;
	closePS=closePS_;
	coutmode=coutmode_;
	myexit=myexit_;
	viderbuff=viderbuff_;
	Commentaire=Commentaire_;
	NoirEtBlanc=NoirEtBlanc_;
	MettreDansPostScript=MettreDansPostScript_;
	getprog=getprog_;
	setgrey=setgrey_;
	getgrey=getgrey_;
  return 1;
}
//  to call the init function before main 
static int init_ff_graph_ptr_func_call = init_ff_graph_ptr_func();

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
#define SetColorTable1 SetColorTable1_
#define GetHeigthFont GetHeigthFont_
//#define compile compile_
//#define compileString compileString_
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
// end wrapping ----	    
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

