// -*- Mode : c++ -*-
//
// SUMMARY  :      
// USAGE    :        
// ORG      : 
// AUTHOR   : Frederic Hecht
// E-MAIL   : hecht@ann.jussieu.fr
//

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

// ALH - javascript output
#ifdef FFJS
#ifndef NATIVEFFJS // [[file:~/ffjs/Makefile::NATIVEFFJS]]
#include <emscripten.h>
#endif // NATIVEFFJS
#endif // FFJS

#define FF_GRAPH_SET_PTR
#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <cassert>
using namespace std;
#include "rgraph.hpp"

#include <sstream> // add by fujiwara
#include "pdf.h" // add by fujiwara

#include "error.hpp"
#ifdef macintoshxx
#include <ConditionalMacros.h>
#include <unix.h>
#else
#include <sys/stat.h>
#endif
#ifdef xTARGET_CARBON
#include <Navigation.h>

int pStrCopy (StringPtr p1, char * p2)
/* copies a pascal string `p1 into a C string */
{
	int len,i;
	
	len = (*p1++) %256;
	for(i=1;i<=len;i++) *p2++=*p1++;
	*p2 = 0;
	return 0;
}

int getprog(char* fn,int argc, char** argvptr)
{
   OSErr anErr;

  NavDialogOptions dialogOptions;
  NavReplyRecord reply;
  
  anErr=NavGetDefaultDialogOptions(&  dialogOptions);
  if( anErr != noErr)  return -1;
  anErr =NavChooseFile(0,&reply,&dialogOptions,0,0,0,0,0) ;  
  if (anErr == noErr && reply.validRecord)
   {
                //  Deal with multiple file selection
                long    count;
                
                anErr = AECountItems(&(reply.selection), &count);
                // Set up index for file list
                if (anErr == noErr)
                {
                    long index;
                    
                    for (index = 1; index <= count; index++)
                    {
                        AEKeyword   theKeyword;
                        DescType    actualType;
                        Size        actualSize;
                        FSSpec      documentFSSpec;
                        
                        // Get a pointer to selected file
                        anErr = AEGetNthPtr(&(reply.selection), index,
                                            typeFSS, &theKeyword,
                                            &actualType,&documentFSSpec,
                                            sizeof(documentFSSpec),
                                            &actualSize);
                        if (anErr == noErr)
                         {  
                           anErr = HSetVol(0,documentFSSpec.vRefNum,documentFSSpec.parID);
                           pStrCopy(documentFSSpec.name, fn);
                         }
                    }
                }
                //  Dispose of NavReplyRecord, resources, descriptors
                anErr = NavDisposeReply(&reply);
     
   }
	return (2);
}
#else
#include "getprog-unix.hpp"
#endif


template<class T> inline T Min (const T &a,const T &b){return a < b ? a : b;}
template<class T> inline T Max (const T &a,const T & b){return a > b ? a : b;}


static  long  cube6[7][3] ={ { 65535,32000,32000},{ 65535, 65535,0},{0, 65535,0},{0, 65535, 65535},{0,0, 65535}
     , { 65535,0, 65535},{ 32000,0,0} }; 
static  long grey6[2][3] ={ {65534,65534,65534},{0,0,0} }; 

static bool grey=false;
static FILE *psfile = 0;
static FILE *psfile_save = 0;
static bool pdffile = false; // add by fujiwara
static char *pdffile_name = nullptr; // add by fujiwara
static float pdf_s = 1; // add by fujiwara
static std::stringstream pdffile_content; // add by fujiwara
static FILE *svgfile = 0; // add by fujiwara
static FILE *svgfile_save = 0; // add by fujiwara
static float svg_s = 1; // add by fujiwara
static int svg_r = 0; // add by fujiwara
static int svg_g = 0; // add by fujiwara
static int svg_b = 0; // add by fujiwara
static int svg_lw = 1; // add by fujiwara
static int LastColor=2;  //  pour est en couleur par defaut

const float fMinPixel = -32000;
const float fMaxPixel = +32000;
#define reel float 
 typedef struct XColor { 
  unsigned short red,green,blue;
} XColor;
static XColor *colortable;
static int ncolortable,fcolor;
static reel echx,echy,rxmin,rxmax,rymin,rymax;
static int  lacouleur=0, width, height, currx=0, curry=0;
#define call(i) i
static int INITGRAPH=0;
void myend();
void myend()
{
 if (INITGRAPH)
   closegraphique();
  cout << "the end" <<endl;
//  ExitToShell();
}
void myexit(int err) { 
  cout << " The End err=" << err << endl;
  exit(err);}

const char * edpfilenamearg=0;	 	
bool  waitatend=false;
bool  consoleatend=false;

#ifdef FREEFEM
#include <fstream.h>
#include <new.h>

void out_of_memory ();
void myexit(int );
void compile(char *fname);


int main (int argc, char **argv)
{
  atexit(myend);

  int OPTION = 0;
  if (argc == 2)
    {
       printf ("PROGRAM  FreeFem 1.0 %s \n",argv[1]);
       initgraphique();
       compile (argv[1]);
       closegraphique();
    }
  else
    printf ("To launch freefem you must type freefem  and a file name\n");
  return 0;
}
#else
void doatexitff();
void doatexitff()
{
#ifdef _WIN32
  bool err=true;
  if(edpfilenamearg && consoleatend)
	{
	cout << " try getConsole " << edpfilenamearg << endl; 
	string fn = edpfilenamearg;
	err=GetConsoleBuff(fn);
        }	
  if( waitatend &&  err)
    {
      char c;  
      cout << "wait enter ? ";
      cin.get();
    }
#endif

}

extern int mymain(int argc,char **argv);

// <<main>>

// ALH - 4/4/15 - main() is not desired if FF is used as a library

#ifndef WITHOUT_MAIN

// ALH - 24/2/15 - the main entry point for javascript is on the HTML page so we need to give it a different name here.
// Built by [[file:~/fflib/Makefile::FFLIB_MAIN]] and called by [[file:~/fflib/fflib.cpp::calling_fflib_main]].

#ifdef FFLIB_MAIN
int fflib_main(int argc,char **argv)
#else
int main (int argc, char **argv)
#endif
{
  atexit(doatexitff);
  int ret=mymain(argc,argv); // <<calling_mymain>>
  return ret;
}

#endif // WITHOUT_MAIN

#endif // FREEFEM

void message(char *s);
void message(char *s)
{  printf("%s	\n",s);}

void erreur(char *s)
{ message(s); exit(0);}
void *safecalloc(size_t nb, size_t  size);
void *safecalloc(size_t nb, size_t  size)
{
  void* p=NULL;
  p = calloc(nb, size);
  if (p == NULL) printf("Run out of Memory!\n");
  return p;
}
void safefree(void** f);
void safefree(void** f)
{
  if(*f){ free((char*) *f); *f=NULL;}
}
void rflush();
void rflush()
{
}


int LaCouleur() {return  lacouleur;}

void couleur(int c)
{ 
  if ( lacouleur == c) // small optim
    return;

  c= c > LastColor ? 1 : c; // c=Min(c,LastColor); pour noir et blanc
 lacouleur = c;
    float r=1,g=1,b=1;
    if (colortable) {
      if (c>0 && c < ncolortable)
	{
	  r =  (float) colortable[c].red /65535.;
	  g =  (float) colortable[c].green /65535.;
	  b =  (float) colortable[c].blue /65535.;
	}
    }
    else if (c!=0)
      r=g=b=0;
    if(psfile)
    fprintf(psfile,"%.3f %.3f %.3f C\n",r,g,b);

    // add by fujiwara
    if( pdffile ){
	pdffile_content << r << ' ' << g << ' ' << b << " RG" << std::endl;
	pdffile_content << r << ' ' << g << ' ' << b << " rg" << std::endl;
    }
    if( svgfile ){
	svg_r = static_cast<int>(r*256);
	svg_g = static_cast<int>(g*256);
	svg_b = static_cast<int>(b*256);
    }

#ifdef FFJS_GRAPH
    // ALH - <<ffjs_couleur>> javascript graph [[file:~/ffjs/main.js::ffjs_couleur]]
    EM_ASM_DOUBLE({ffjs_couleur($0,$1,$2)},r,g,b);
#endif
}

extern void DefColor(float & r, float & g, float & b,
              int k,int nb, bool hsv,bool ggrey,int nbcolors,float *colors);

static XColor DefColorSansG( int k,int nb, bool hsv,bool ggrey,int nbcolors,float *colors)
{
 XColor C;
 float r,g,b;
 DefColor(r,g,b,   k,nb,hsv,ggrey,nbcolors,colors);
 C.red= (short unsigned int) (65535*r);
 C.green=(short unsigned int)(65535*g);
 C.blue= (short unsigned int) (65535*b);
 return C;
} 
void SetColorTable1(int nb,bool hsv,int nbcolors,float *colors)
{
  static bool greyo = !grey;
  static float *colorso =0;
  if(!INITGRAPH) return;
   if (ncolortable == nb && greyo == grey && colorso == colors ) return;// optim
   greyo = grey;
   colorso=colors;
   if (nbcolors && nb>2) 
     { 
       if(colortable) delete [] colortable;
       colortable = new XColor[nb];
       ncolortable = nb;
       if(LastColor>1) LastColor=nb-1;
        for (int i0=0;i0<nb;i0++)
         {  
           colortable[i0]=DefColorSansG(i0,nb,hsv,grey,nbcolors,colors);           
          }
        
       }
     else 
      ncolortable  =0;
}
void SetColorTable(int nb)
{
   if (fcolor  && nb>2 && nb < 256) 
     { 
       nb = Max(nb,8);
       if (ncolortable == nb) return;// optim
       if(colortable) delete [] colortable;
       colortable = new XColor[nb];
       ncolortable = nb;
       if(LastColor>1) LastColor=nb-1;
      // cout << "SetColorTable "<< nb << endl;
       int k=0;
       colortable[k].red= 65535;
       colortable[k].green= 65535;
       colortable[k].blue= 65535;
       k++;
       colortable[k].red=0;
       colortable[k].green=0;
       colortable[k].blue=0;
       k++;
       nb = nb -2;
       for (long i0=0;i0<nb;i0++,k++)
         {  
	   //     long  i1 = nb - i0;
	   long  i6 = i0*6;
	   long  j0 = i6/nb;// in 0..6
	   long  j1 = j0+1;// in 1..6
	   long  k0 = i0 - (nb*j0)/6L;
	   long  k1 = (nb*j1)/6L-i0;
	   long  kk = k0+k1;
	   //	     cout <<k << " " << i0 << " " << j0 << " " << j1 << " " << k0 << " " << k1  << endl;
           if ( kk <= 0)
	     { cerr << kk << " " << nb << " " << k0 << " " << k1 << " " << endl;
	     assert(kk);
             }
      /*     colortable[k].red   = (unsigned short) ((long) (cube6[j1][0]*k0+cube6[j0][0]*k1)/kk);
           colortable[k].green = (unsigned short) ((long)  (cube6[j1][1]*k0+cube6[j0][1]*k1)/kk);
           colortable[k].blue  = (unsigned short) ((long) (cube6[j1][2]*k0+cube6[j0][2]*k1)/kk);*/
	   if (! grey)
         {
           colortable[k].red   = (unsigned short) ((long) (cube6[j1][0]*k0+cube6[j0][0]*k1)/kk);
           colortable[k].green = (unsigned short) ((long)  (cube6[j1][1]*k0+cube6[j0][1]*k1)/kk);
           colortable[k].blue  = (unsigned short) ((long) (cube6[j1][2]*k0+cube6[j0][2]*k1)/kk);
	   }
          else 
           {
           kk=nb-1;
           k1 =  i0;
           k0 = nb - i0 -1;
           j0=1;
           j1=0;
           colortable[k].red   = (unsigned short) ((long) (grey6[j1][0]*k0+grey6[j0][0]*k1)/kk);
           colortable[k].green = (unsigned short) ((long)  (grey6[j1][1]*k0+grey6[j0][1]*k1)/kk);
           colortable[k].blue  = (unsigned short) ((long) (grey6[j1][2]*k0+grey6[j0][2]*k1)/kk);
           }
           
           
           assert(k<ncolortable);
	   //   cout <<colortable[k].pixel 
	   //	<< " r=" <<  colortable[k].red 
	   //	<< " g=" <<  colortable[k].green
	   //	<< " b=" <<colortable[k].blue <<endl;
          }
  }


// a faire 
}
void FlushEvent()
{
} 

void initgraphique()
{
  colortable=0;
  ncolortable=0;
  LastColor=2;// En couleur par default
  fcolor=1;  /* des couleurs */
  SetColorTable(8);
  INITGRAPH = 1;
  width = 10000;// change FH mai 2012 to have more precis graphic for F.  Ortegon 
  height =  7071; // <<aspect_ratio>>  \sqrt(2)
}

void closegraphique()
{
  if (INITGRAPH)
    {
      INITGRAPH = 0;
      delete [] colortable;
      closePS();
    }
}

void cadre(reel xmin,reel xmax,reel ymin,reel ymax)
{
  rxmin = xmin;
  rxmax = xmax;
  rymin = ymin;
  rymax = ymax;

  echx = width / (xmax - xmin);
  echy = height / (ymax - ymin);
}

void getcadre(reel &xmin,reel &xmax,reel &ymin,reel &ymax)
{
  xmin = rxmin;
  xmax = rxmax;
  ymin = rymin;
  ymax = rymax;

}


int InRecScreen(reel x1, reel y1,reel x2, reel y2)
{  

  return (Max(x1,x2)>= rxmin) && (Min(x1,x2) <= rxmax) && (Max(y1,y2) >= rymin) && (Min(y1,y2)  <= rymax);
}
int InPtScreen( reel x, reel y)
{
  return (x >= rxmin) && (x <= rxmax) && (y >= rymin) && (y <= rymax);
}


float scali(int i);
float scalj(int i);
float scalx(int i);
float scaly(int i);
float scali(int i)
{
  return i/echx  + rxmin;
}
float scalj(int j)
{
  return -j/echy  + rymax;
}
int scalx(reel x)
{
  return (int) Min(fMaxPixel,Max(fMinPixel,((x - rxmin) * echx)));
} 
int scaly(reel y)
{
  return (int)Min(fMaxPixel,Max(fMinPixel,((rymax - y) * echy)));
}

void pointe(reel , reel )
{
}

// <<rmoveto>>
void rmoveto(reel x, reel y)
{
  currx = scalx(x);
  curry = scaly(y);

#ifdef FFJS_GRAPH
  // ALH - <<ffjs_rmoveto>> javascript graph [[file:~/ffjs/ffapi.js::ffjs_rmoveto]]
  EM_ASM_INT({ffjs_rmoveto($0,$1)},currx,curry);
#endif
}

void rlineto(reel x, reel y)
{
  int newx = scalx(x), newy = scaly(y);
  if (psfile)
    fprintf(psfile,"%d %d %d %d L\n",currx, height-curry, newx, height-newy);
  // add by fujiwara
  if ( pdffile ){
      pdffile_content << currx*pdf_s << ' ' << (height-curry)*pdf_s << " m "
      		      << newx*pdf_s << ' ' << (height-newy)*pdf_s << " l S" << std::endl;
  }
  if ( svgfile ){
    fprintf(svgfile,"<line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" stroke=\"rgb(%d,%d,%d)\" stroke-width=\"%d\" />\n",
	    currx*svg_s, curry*svg_s, newx*svg_s, newy*svg_s, svg_r, svg_g, svg_b, svg_lw);
  }
  // add by fujiwara
  currx = newx; curry = newy;

#ifdef FFJS_GRAPH
  // ALH - <<ffjs_rlineto>> javascript graph [[file:~/ffjs/main.js::ffjs_rlineto]]
  EM_ASM_INT({ffjs_rlineto($0,$1)},newx,newy);
#endif
}

void cadreortho(reel centrex, reel centrey, reel rayon)
{

  if (height < width)
  {
    rymin = centrey - rayon;
    rymax = centrey + rayon;
    echx = echy= height / (2 * rayon);
    rxmin= centrex - width / (2 * echx);
    rxmax= centrex + width / (2 * echx);
  }
  else
  {
    rxmin = centrex - rayon;
    rxmax = centrex + rayon;
    echx = echy = width / (2 * rayon);
    rymin = centrey - height / (2 * echy);
    rymax = centrey + height / (2 * echy);
  }
}

void plotstring (const char *  string)
{ //int l = strlen(string);
  if(psfile) fprintf(psfile,"(%s) %d %d  S\n",string,currx,height-curry);
 
  // add by fujiwara
  if( pdffile ){
      pdffile_content << "BT /F1 " << 9 << " Tf" << std::endl; // 9*pdf_s : text font size
      pdffile_content << "1 0 0 1 " << currx*pdf_s << ' ' << (height-curry)*pdf_s << " Tm" << std::endl;
      pdffile_content << "(" << string << ") Tj ET" << std::endl;
  }
  if( svgfile ){
      fprintf(svgfile,"<text x=\"%.2f\" y=\"%.2f\">%s</text>\n",currx*svg_s,curry*svg_s,string);
  }
  // add by fujiwara
 
#ifdef FFJS_GRAPH
  // ALH - <<ffjs_plotstring>> javascript graph - Send the string character by character to
  // [[file:~/ffjs/main.js::ffjs_plotstring]] because there is no string parameter at the moment
  // [[http://kripken.github.io/emscripten-site/docs/api_reference/emscripten.h.html#c.EM_ASM_]]
  
  for(const char *c=string;*c;c++)EM_ASM_INT({ffjs_plotstring($0)},*c);
  EM_ASM(ffjs_plotstring(0));
#endif
}

void showgraphic()
{
}

void penthickness(int pepais)
{
  if (psfile) fprintf(psfile,"%d setlinewidth\n",pepais*2);
  // add by fujiwara
  if ( pdffile ){
      pdffile_content << pepais << " w" << std::endl;
  }
  if ( svgfile ){
      svg_lw = pepais;
  }
  // add bu fujiwara

#ifdef FFJS_GRAPH
  // ALH - <<ffjs_penthickness>> javascript graph [[file:~/ffjs/main.js::ffjs_penthickness]]
  EM_ASM_INT({ffjs_penthickness($0)},pepais);
#endif
}

void x11linsrn(int * ,int * ,int * ,int * );
void x11linsrn(int * ,int * ,int * ,int * )
  //int *x1,*x2,*y1,*y2;
{   
}

   
void viderbuff()
{
}


void cercle(reel , reel , reel );
void cercle(reel , reel , reel )
{
  //int r = (int) (rayon * echx);
}
void reffecran(){
#ifdef FFJS_GRAPH
  // ALH - <<ffjs_graphstart>> javascript graph [[file:~/ffjs/main.js::ffjs_graphstart]]
  EM_ASM(ffjs_graphstart());
#endif
}

void fillpoly(int n, float *poly)
{
  int i;
  if (psfile) 
    {
      fprintf(psfile,"bF ");
      for (i=0;i<n;i++)
	fprintf(psfile,"%d %d ", scalx(poly[2*i]),height-scaly( poly[2*i+1]));
      fprintf(psfile,"eF\n");
    }
  // add by fujiwara
  if ( pdffile ) {

      i=0;
      pdffile_content << scalx(poly[2*i])*pdf_s << ' ' << (height-scaly( poly[2*i+1]))*pdf_s << " m ";
      for (i=1;i<n;i++)
	  pdffile_content << scalx(poly[2*i])*pdf_s << ' ' << (height-scaly( poly[2*i+1]))*pdf_s << " l ";
      pdffile_content << "f" << std::endl;
  }
  if ( svgfile ) 
    {
      fprintf(svgfile, "<polygon points=\"");
      for (i=0;i<n;i++)
	fprintf(svgfile,"%.2f,%.2f ", scalx(poly[2*i])*svg_s,scaly( poly[2*i+1])*svg_s);
      fprintf(svgfile,"\" stroke=\"rgb(%d,%d,%d)\" stroke-width=\"%d\" fill=\"rgb(%d,%d,%d)\" />\n",
	      svg_r,svg_g,svg_b, svg_lw, svg_r,svg_g,svg_b);
    }
  // add by fujiwara
   
#ifdef FFJS_GRAPH
  // ALH - <<ffjs_fillpoly>> javascript graph [[file:~/ffjs/main.js::ffjs_fillpoly]]
  EM_ASM_INT({ffjs_fillpoly_begin($0,$1)},scalx(poly[0]),scaly(poly[1]));
  for(i=1;i<n;i++)EM_ASM_INT({ffjs_fillpoly_next($0,$1)},scalx(poly[2*i]),scaly(poly[2*i+1]));
  EM_ASM(ffjs_fillpoly_close());
#endif
}

//----------------------------------------------------------------------
// add by fujiwara
//----------------------------------------------------------------------
void openPDF(const char *filename )
{
  if(pdffile) closePDF();

  if( pdffile_name != nullptr ){
      delete [] pdffile_name;
      pdffile_name = nullptr;
  }

  pdffile_name = new char [ strlen(filename)+1 ];
  strcpy( pdffile_name, filename );

  pdffile_content.str(""); // clear
  pdffile_content.clear( std::stringstream::goodbit );

  pdffile_content.setf( std::ios::fixed );
  pdffile_content.precision( 3 );

  const int widthA4PDF = 596;
  pdf_s = static_cast<float>(widthA4PDF) / width;

  pdffile = true;
  pdffile_content << "q" << std::endl; // gsave
  
  return;
}
void closePDF(void)
{
  if(pdffile) {

      std::string PDFTitle = "plot() by FreeFem++";
      std::string AppName = "FreeFem++ v" + StrVersionNumber();

      SimplePDF_FF pdf( pdffile_name, PDFTitle.c_str(), AppName.c_str() );

      const int widthPDF  = static_cast<int>( width * pdf_s );
      const int heightPDF = static_cast<int>( height * pdf_s );

      pdffile_content << "Q" << std::endl;

      pdf.addPage( pdffile_content, widthPDF, heightPDF );

      pdffile = false;
  }

  if( pdffile_name != nullptr ){
      delete [] pdffile_name;
      pdffile_name = nullptr;
  }

  return;
}
void openSVG(const char *filename )
{
  if(svgfile_save) closeSVG();

  const int  widthA4PS = 596;
  //const int heightA4PS = 842;
  svg_s = static_cast<double>(widthA4PS)/width;

  char ffff[32];
  int count = 0;
  if(!filename){
    bool notfound;
    do {
      struct stat buf;
      snprintf(ffff,32,"rgraph_%.3d.svg",count++);
      volatile int r = stat(ffff,&buf) ;
      notfound = (r != 0);
      if( count > 1000 ) break;
    } while ( !notfound );
  }   

  const char *fsvg (filename?filename:ffff);

  svgfile=fopen(fsvg,"w");

  if(svgfile) {
    svgfile_save=svgfile;
    fprintf(svgfile,"<?xml version=\"1.0\" standalone=\"yes\"?>\n\n");
    fprintf(svgfile,"<!-- Creator: FreeFem++ v%s -->\n\n", StrVersionNumber().c_str());
    fprintf(svgfile,"<svg xmlns=\"http://www.w3.org/2000/svg\"\n");
    fprintf(svgfile,"\txmlns:xlink=\"http://www.w3.org/1999/xlink\"\n");
    fprintf(svgfile,"\twidth=\"%dpx\" height=\"%dpx\">\n",
	    static_cast<int>(width*svg_s), static_cast<int>(height*svg_s));
  }
  else
  {
    cerr << " Err opening SVG file " << fsvg << endl;
  }
  return;
}
void closeSVG(void)
{
  if(svgfile_save) {
    fprintf(svgfile_save,"</svg>\n");
    fclose(svgfile_save);
  }
  svgfile_save=0;
  svgfile=0;

  return;
}
//----------------------------------------------------------------------
// add by fujiwara end
//----------------------------------------------------------------------


int  execute (const char * str)
{ 
 if(verbosity)
     cout << "exec: " << str << endl;
 return  system(str);
}

char Getijc(int *x1,int *y1);
char Getijc(int *x1,int *y1)
{
  //char char1;
  *x1=0;
  *y1=0;
  //cout << "entre un caractere ? ";
  //cin >>char1 ;
  return 0;// char1;
}
   
char Getxyc(float &x,float &y)
{ 
  //  cout << " in Getxyc" << endl;
  char c;
  int i,j;
  c = Getijc( &i,&j);
  x = scali(i);
  y = scalj(j);
  //  cout << " out  Getxyc" << x << " " << y << " " << c << endl;

  return c;
}

/// <<rattente>>
void rattente(int )
{
}
 void GetScreenSize(int &ix,int &iy)
{
  ix = width;
  iy = height;
}

// <<openPS>>

void openPS(const char *filename )
{ 
  char ffff[32];
  int count=0;
  if(psfile_save) closePS();
  time_t t_loc;
  int  widthA4PS=596;
  //int heightA4PS=842;
  float s= (double)widthA4PS/width;
  char  username[10];
  /*if (!cuserid(username)) */ strcpy(username,"inconnue");
  time(&t_loc);
  bool notfound;
  if( !filename) 
   do {
      struct stat buf;
      snprintf(ffff,32,"rgraph_%.3d.ps",count++);
      volatile int r= stat(ffff,&buf) ;
      notfound = r !=0;
      if(count>1000) break;
    } while ( !notfound );
   

  const char *fps (filename?filename:ffff);


  psfile=fopen(fps,"w");

#ifdef FFJS_GRAPH

  // <<listgraphs>> ALH - extract the list of graphs generated by FF during a run
  // [[file:~/ffjs/ffapi.js::ffjs_listgraphs]]. Send the file name character by character to
  // [[file:~/ffjs/ffapi.js::ffjs_listgraphs]] because there is no string parameter at the moment
  // [[http://kripken.github.io/emscripten-site/docs/api_reference/emscripten.h.html#c.EM_ASM_]]
  
  for(const char *c=fps;*c;c++)EM_ASM_INT({ffjs_listgraphs($0)},*c);
  EM_ASM(ffjs_listgraphs(0));
#endif

  if(psfile) {
    psfile_save=psfile;
    fprintf(psfile,"%%!PS-Adobe-2.0 EPSF-2.0\n%%%%Creator: %s\n%%%%Title: FreeFem++\n","user");
    fprintf(psfile,"%%%%CreationDate: %s",ctime(&t_loc));
    fprintf(psfile,"%%%%Pages: 1\n");
    fprintf(psfile,"%%%%BoundingBox:       0 0 %d %d\n",int(width*s),int(height*s));
    fprintf(psfile,"%%%%EndComments\n");
    fprintf(psfile," /L {newpath moveto lineto stroke} def\n");
    fprintf(psfile," /C {setrgbcolor} def\n");
    fprintf(psfile," /rec {newpath 4 copy 8 1 roll moveto 3 -1 roll lineto 4 2 roll exch lineto lineto closepath} def\n");
    fprintf(psfile," %f %f  scale \n",s,s);
    fprintf(psfile," 0 %d 0 %d rec clip\n",int(width),int(height));
    fprintf(psfile," /Helvetica findfont %d scalefont setfont\n",int(9/s));
    fprintf(psfile," /S {moveto show} def\n");
    fprintf(psfile," /bF  { mark} def \n");
    fprintf(psfile," /eF {newpath moveto counttomark 2 idiv {lineto} repeat closepath fill cleartomark} def\n");
    fprintf(psfile," /P { /yy exch def /xx exch def   xx xx 1 add yy yy 1 add  rec  fill } def\n");
    fprintf(psfile," 2 setlinewidth\n");
  }
  else 
    cerr << " Err opening postscript file " << fps << endl;
}
void closePS(void)
{
  if(psfile_save) {
    fprintf(psfile_save,"showpage\n");
    fclose(psfile_save);
  }
  psfile_save=0;
  psfile=0;

#ifdef FFJS_GRAPH
  // ALH - <<ffjs_graphdone>> javascript graph [[file:~/ffjs/main.js::ffjs_graphdone]]
  EM_ASM({ffjs_graphdone()});
#endif
}
 void coutmode(short )  {}
// bof bof --- 
 float  GetHeigthFont()
{ 
     double  widthA4PS=596;
    float s=widthA4PS/width; 
  return 5.5/s/echy;
}
  void Commentaire(const char * c)  
  {
  if(psfile)   {
    fprintf(psfile,"%% %s\n",c);
   }
  // add by fujiwara
  if( pdffile ) {
      //fprintf(pdffile,"%% %s\n",c);
   }
  if( svgfile ) {
    fprintf(svgfile,"%% %s\n",c);
   }
  // add by fujiwara
  }
  void NoirEtBlanc(int NB)
  {
    if(NB) LastColor=1;
    else LastColor=ncolortable?ncolortable:2;
  }
 
  void MettreDansPostScript(int in)
   {
     if(in)  psfile=psfile_save;     
     else   psfile=0;
   }
  // add by fujiwara
   void MettreDansPDF(int in) // put into PDF
   {
     if(in)  pdffile=true;
     else   pdffile=false;
   }
   void MettreDansSVG(int in) // put into SVG
   {
     if(in)  svgfile=svgfile_save;
     else   svgfile=0;
   }
   // add by fujiwara

static void     FillRect(float x0,float y0, float x1, float y1)
 {
     float r[8];
     r[0]=x0;r[1]=y0;
     r[2]=x1;r[3]=y0;
     r[4]=x1;r[5]=y1;
     r[6]=x0;r[7]=y1;
     fillpoly(4,r);
 }
int PutLevel(int lineno, float xf, int col);
int PutLevel(int lineno, float xf, int col)
{
  float xmin,xmax,ymin,ymax;
  getcadre(xmin,xmax,ymin,ymax);
  float xleft = xmax - (xmax-xmin)*0.1;
  float ytop  = ymax;
  float ydelta = (ymax-ymin)/40;
  ydelta=GetHeigthFont();
  xleft = xmax - 6*ydelta;  
  ytop -= ydelta*(col+2);
  couleur(col);
  FillRect(xleft+ydelta/8.,ytop+ydelta/8.,xleft+ydelta*7./8.,ytop+ydelta*7./8.);
  rmoveto(xleft+ydelta*1.4,ytop+ydelta/4);
  char buf[30];
  snprintf(buf,30,"%g",xf);
  couleur(1);
  plotstring(buf);

   return lineno;
}
 void ShowHelp(const char * s,int k)
{
  if(k) {
    MettreDansPostScript(0);
    MettreDansPDF(0); // add by fujiwara
    MettreDansSVG(0); // add by fujiwara
    couleur(1);
    float xmin,xmax,ymin,ymax;
    getcadre(xmin,xmax,ymin,ymax);
    rmoveto(xmin+(xmax-xmin)/100,ymax-(k)*(ymax-ymin)/30);
    plotstring(s);
    MettreDansPostScript(1);
    MettreDansPDF(1); // add by fujiwara
    MettreDansSVG(1); // add by fujiwara
       //  couleur(1);	
  }
}

  void setgrey(bool gg ){grey=gg;}
  int getgrey(){ return grey;}

class Grid;
void SaveMesh(Grid &);
void SavePlot(int , Grid& , double *);
void SavePlot(int , Grid& , float *);

void SaveMesh(Grid &){}
void SavePlot(int , Grid& , double *){}
void SavePlot(int , Grid& , float *){}

