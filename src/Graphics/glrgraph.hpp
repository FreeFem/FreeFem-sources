/**************DO NOR REMOVE THIS BANNER***************/
/*  FreeFEM : Language for a Finite Element Method    */
/*  -------    Release 1.0:  June 1994.               */
/*  Authors: D. Bernardi, Y. Darmaillac F. Hecht,     */
/*           O. Pironneau                             */
/*  You may copy freely these files and use it for    */
/* teaching or research. These or part of these may   */
/* not be sold or used for a commercial purpose with- */
/* out our consent : fax (33)1 44 27 44 11            */
/* (fax)    Olivier.Pironneau@ann.jussieu.fr          */
/******************************************************/
//#define FREEFEM
//  AGL  apple 
//  XGL   X11  
//  WGL   window (a faire) 
#ifdef AGL
#define TARGET_API_MAC_CARBON 1
#define CALL_IN_SPOCKETS_BUT_NOT_IN_CARBON 1
#include <Carbon/Carbon.h>
    
#include <AGL/agl.h>
#endif
#ifdef XGL
#include <GL/glx.h>
#include <X11/cursorfont.h>
#include <X11/keysymdef.h>
#endif

#include <GL/gl.h>
#include <GL/glu.h>

#include <sys/stat.h>

#include "error.hpp"
#include <fstream>
#include <iostream>
#include <stdlib.h>
using namespace std;

#include <errno.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "rgraph.hpp"
#include <time.h>
#include <setjmp.h>
#include <time.h>

int currx=0,curry=0;
static FILE *psfile = 0;
static FILE *psfile_save = 0;

#ifdef AGL
static	AGLPixelFormat fmt;
static	AGLContext ctx;

int pStrCopy (StringPtr p1, char * p2);
#endif
#ifdef XGL
static  Display *dpy;
static  Window win;
static  XSizeHints size_hints;
//static  GC gc;
static  XFontStruct *font_info;
GLXContext cx;
int stensize;
static Cursor cursor_watch,cursor_arrow;
static int shift, control,shiftlock,alt;
static GLuint basefont; 
#endif


extern long verbosity;  // level off printing


#ifdef FREEFEM
void myenviron ()
{
  cout << "FreeFEM error: operator new failed; not enough memory" << endl;
  if (myenviron)
   longjmp(myenvironj,1);
  exit(2);
}
//  pour imprimer la version   FH 
#define STRING(i) #i
#include <new.h>

jmp_buf myenvironj;
static int  myenviron = 0;

void out_of_memory ();
void NEW_HANDLER (void);
void compile(char *fname);
float scali(int i);
float scalj(int j);
void execute(char* what);
char Getijc(int & x,int & y);
int DoMouseDown (int windowPart, WindowPtr whichWindow, EventRecord *myEvent);
 void NEW_HANDLER (void){  set_new_handler (&myenvironj);}
#endif

static int nbcolor;
static int ncolortable;
static int LastColor; // LastColor=1 => Noir et Blanc 

#ifdef AGL
#define	ours(w)		(w==grafWindow0)
static WindowPtr	 grafWindow0;
static GrafPtr          grafPort0;
static	Rect		boundsRect;
static CursHandle  CrossCurseur ;
static CursHandle  WatchCurseur ;
static  Pattern  white,black;
#else
struct RGBColor {
  unsigned short      red;                    /*magnitude of red component*/
  unsigned short      green;                  /*magnitude of green component*/
  unsigned short      blue;                   /*magnitude of blue component*/
};
#endif

template<class T> inline T Min (const T &a,const T &b){return a < b ? a : b;}
template<class T> inline T Max (const T &a,const T & b){return a > b ? a : b;}

static bool grey=false;
static  int cube6[7][3] ={ {65534,0,0},{65534,65534,0},{0,65534,0},{0,65534,65534},{0,0,65534}
     , {65534,0,65534},{65534,0,0} }; 
static  int grey6[2][3] ={ {65534,65534,65534},{0,0,0} }; 

char errbuf[255];
static int INITGRAPH=0;
static float aspx, aspy, echx,echy,ech,rxmin,rxmax,rymin,rymax;
static int carre, lacouleur;
static	GLuint fontList;



static int width,height;

static RGBColor * colortable;
int getcolor();
void putpixel(int ix,int iy, int couleur);
int scalx(float x);
int scaly(float y);
void thisexit();


void DrawCStringGL (const char * cstrOut, GLuint fontList)
{
	GLint i = 0;
	glRasterPos3d(currx,height-curry,0);
	while (cstrOut [i])
		glCallList (fontList + cstrOut[i++]);
}


#ifdef AGL
void InitMac();
// --------------------------------------------------------------------------------------------------------------

// --------------------------------------------------------------------------------------------------------------
// APPLE EVENT SUPPORT ROUTINES
// --------------------------------------------------------------------------------------------------------------

int pStrCopy (StringPtr p1, char * p2)
/* copies a pascal string `p1 into a C string */
{
	int len,i;
	
	len = (*p1++) %256;
	for(i=1;i<=len;i++) *p2++=*p1++;
	*p2 = 0;
	return 0;
}

void InitMac()
{
	BitMap	screenBitMap;
	Rect	screenBits;
	Cursor theArrow;
	GetQDGlobalsScreenBits(&screenBitMap);
	screenBits = screenBitMap.bounds;
	SetCursor(GetQDGlobalsArrow(&theArrow));

}
class InitilisationMac {
  static int init;
  public:
   InitilisationMac(){ InitMac();}
};

static InitilisationMac Initmac; // to call InitMac

int getprog(char* fn,int  argc, char** argv)
{ 
 if (argc > 2) 
  {
    for (int i=1; i<argc-1;i++)
     if  (strcmp(argv[i],"-f")==0 )
       {
          strcpy(fn,argv[i+1]);
          printf(" file : %s\n",fn);
          return 1;
          
       }
      else if(strcmp(argv[i],"-v")==0) 
       {  
         verbosity=atoi(argv[++i]);
       }
    printf(" if more than 1 argument then  use -f filename");
    return 0; 
  }
 else if (argc > 1) 
  { 
    strcpy(fn,argv[1]);
    printf(" file : %s\n",fn);
  }
 else 
  {   OSErr anErr;


  FSRef  fsRef;
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
                                            typeFSRef, &theKeyword,
                                            &actualType,&fsRef,
                                            sizeof(fsRef),
                                            &actualSize);
                        if (anErr == noErr)
                         {  
                          anErr=FSRefMakePath(&fsRef,(UInt8*)fn,256);
                          if ( anErr == noErr )
                           {
                           cout <<  "Path : " << fn << endl;
                           char * ff=fn,*fff=0;
                           while ( *ff)
                            { 
                              if (*ff=='/') fff=ff;
                              ff++;
                            }
                           if (fff) { 
                             *fff=0;                              
                             cout << "chdir to "<< fn << endl;
                            chdir(fn);
                             *fff='/';}
                           }
                           else cout << "Err: "<< anErr << endl;
                           /*
                        anErr = AEGetNthPtr(&(reply.selection), index,
                                            typeFSS, &theKeyword,
                                            &actualType,&documentFSSpec,
                                            sizeof(documentFSSpec),
                                            &actualSize);
                           anErr = HSetVol(0,documentFSSpec.vRefNum,documentFSSpec.parID);
                           pStrCopy(documentFSSpec.name, fn);*/
                         }
                    }
                }
                //  Dispose of NavReplyRecord, resources, descriptors
                anErr = NavDisposeReply(&reply);
     
   }
   else return 0; // erreur cancel
	return (2);
	}
}


//-----------------------------------------------------------------------------------------------------------------------

GLuint BuildFontGL (AGLContext ctx, GLint fontID, Style face, GLint size)
{
	GLuint listBase = glGenLists (256);
	if (aglUseFont (ctx, fontID , face, size, 0, 256, (long) listBase))
	{
		glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);
		return listBase;
	}
	else
	{
		glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);
		glDeleteLists (listBase, 256);
		return 0;
	}
}
#else
int getprog(char* fn,int argc, char **argv)
{
 if (argc > 2) 
  {
    for (int i=1; i<argc-1;i++)
     if  (strcmp(argv[i],"-f")==0)
       {
          strcpy(fn,argv[i+1]);
          printf(" file : %s\n",fn);
          return 1;
       }
    printf(" if more than 1 argument then  use -f filename");
    return 0; 
  }
 else if (argc > 1) 
  { 
    strcpy(fn,argv[1]);
    printf(" file : %s\n",fn);
  }
 else {
   cout << " no argument " << endl;
   return 0;
 }
 return argc;
}

#endif

#ifdef XGL
void  MyXSelectInput(Display * dpy,Window w,int  mask)

{
  XSetWindowAttributes attributes;
  attributes.event_mask = mask;
  XChangeWindowAttributes(dpy, w, CWEventMask, &attributes);
}
#endif

//-----------------------------------------------------------------------------------------------------------------------

void DeleteFontGL (GLuint fontList)
{
	if (fontList)
		glDeleteLists (fontList, 256);
}

#ifdef FREEFEM

void coutmode(short i) 
{ 
   cout <<  flush;
   cerr <<  flush;
;}

void myexit(int err)
{
 if (INITGRAPH)
  {
    rattente(0);
    closegraphique();
  }
 if (err !=0)
    cout << "Error: freefem+ has end with error code " <<err<<endl;
// else cout << "Normal exit 0" << endl;
  if (myenviron)
   longjmp(myenvironj,1);
}
void thisexit(){ myexit();}

int main (int argc, char **argv)
{
  char       *prog;
  char       fname[256];
  argc = getprog (fname, argc, argv);
  atexit(thisexit);
  NEW_HANDLER (); // see dependent system files ({pc,x,mac}rgraph.{h,cpp})
  

  int OPTION = 0;
  if (argc == 2)
    {
        initgraphique();
        if(0==setjmp(myenvironj))
         {  myenviron=1;
	   compile (fname);
	  // cout << "No Error" << endl;
	 }
	myenviron = 0;			
    }
  else
    printf ("To launch freefem you must type freefem  and a file name\n");
  
  return 0;
}
#else
extern int mymain(int argc,char **argv);
const char * StrVersionNumber();

int main (int argc, char **argv)
{
    char * wn = new  char [256];
   for (int i=0;i<256;i++)
     wn[i] = 0;
   strcpy(wn,"   -- FreeFem++ ");
   strcat(wn,StrVersionNumber());
   
   int ret=15;  
   try {                  
          ret=mymain(argc,argv);}
   catch( Error & err) {
                          cerr  << err.what() << endl;                        
                        }
   catch( ...) { cerr << "catch exception ???";}
                        

 return ret;
}
#endif

void message(char *s)
{ 
   printf("%s	\n", s);
}

void erreur(char *s)
{
    cout  << endl;
    cerr << "##Fatal error  :" << s << endl << "exit(1)" <<  endl;
    exit(1);
}

void *safecalloc(long nb, long size)
{
  void* p=NULL;
  p = calloc(nb, size);
  if (p == NULL) 
     erreur("Out of Memory!\n");
  return p;
}

void safefree(void** f)
{
  if(*f)
  { 
    free(*f); 
    *f=NULL;
  }
}

void initgraphique(void)
{
    if(INITGRAPH) return;
 //   cout <<"Initgraphique \n" ;
    fontList=0;
#ifdef AGL 
	BitMap	screenBitMap;
	Rect	screenBits;
	Cursor theArrow;
	GetQDGlobalsScreenBits(&screenBitMap);
	screenBits = screenBitMap.bounds;
	SetCursor(GetQDGlobalsArrow(&theArrow));

	boundsRect.top = 45;
	boundsRect.left = (short) (15 +  (0.35 * screenBits.right));
	boundsRect.bottom = screenBits.bottom -  25;
	boundsRect.right =  screenBits.right-  25;
	if((boundsRect.bottom - boundsRect.top) < (boundsRect.right - boundsRect.left))
		boundsRect.right = boundsRect.left + boundsRect.bottom - boundsRect.top;
	else
		boundsRect.bottom = boundsRect.top + boundsRect.right - boundsRect.left;
	grafWindow0=NewCWindow(0, &boundsRect, "\pFreeFem Graphics",true, 8, (WindowPtr) -1L, true, 0);
	
	//ShowWindow(grafWindow0);
	BringToFront(grafWindow0);
	//SelectWindow(grafWindow0);
	SetPortWindowPort(grafWindow0);
	GetPort(&grafPort0);
	
	height = boundsRect.bottom - boundsRect.top - 10;
	width = boundsRect.right - boundsRect.left -10;
	aspx = boundsRect.right - boundsRect.left -10;
	aspy = boundsRect.bottom - boundsRect.top - 10;

	
	
	GLint attrib[] = { AGL_RGBA, AGL_DOUBLEBUFFER, AGL_NONE };
	
	fmt = aglChoosePixelFormat(NULL, 0, attrib); /* Choose pixel format */

	ctx = aglCreateContext(fmt, NULL); 	/* Create an AGL context */

	aglDestroyPixelFormat(fmt); // pixel format is no longer needed

	aglSetDrawable(ctx, GetWindowPort (grafWindow0)); /* Attach the context to the window */

	{
		EventRecord event;
		WaitNextEvent (everyEvent, &event, 1, NULL);
	}

	aglSetCurrentContext(ctx);
    short int fNum;
   // cout <<" GetFNum \n";
   
	GetFNum("\pGeneva", &fNum);									// build font
	fontList = BuildFontGL (ctx, fNum, normal, 9); 
#endif
#ifdef XGL
 {
    XVisualInfo* vi;
    Colormap cmap;
    XSetWindowAttributes swa;
    static int attrib[] = { GLX_RGBA,
			    GLX_DOUBLEBUFFER,
			    GLX_RED_SIZE, 1,
			    GLX_GREEN_SIZE, 1,
			    GLX_BLUE_SIZE, 1,
			    GLX_DEPTH_SIZE, 16,
			    GLX_STENCIL_SIZE, 4,
			    None };

    /* get a connection */
      dpy = XOpenDisplay(0);
  if (!dpy) 
    {
      cerr << " Erreur openning  dpy " << endl;
      exit(2);
    }

    /* get an appropriate visual */
    vi = glXChooseVisual(dpy, DefaultScreen(dpy), attrib);
    if (vi == NULL) {
	fprintf(stderr, "Can't find a satisfactory visual.  Abort.\n");
	exit(1);
    }
    glXGetConfig(dpy, vi, GLX_STENCIL_SIZE, &stensize);

    /* create a GLX context */
    cx = glXCreateContext(dpy, vi, 0, GL_TRUE);

    /* create a color map */
    cmap = XCreateColormap(dpy, RootWindow(dpy, vi->screen),
                           vi->visual, AllocNone);

    /* create a window */
    swa.colormap = cmap;
    swa.border_pixel = 0;
    swa.event_mask = StructureNotifyMask;
    height = 512;
    width = 512;
	aspx = width;
	aspy =height;
    win = XCreateWindow(dpy, RootWindow(dpy, vi->screen), 0, 0, width, height,
                        0, vi->depth, InputOutput, vi->visual,
                        CWBorderPixel|CWColormap|CWEventMask, &swa);
    XMapWindow(dpy, win);
	glXMakeCurrent(dpy, win, cx); 
    cursor_arrow = XCreateFontCursor(dpy,XC_arrow);
    cursor_watch = XCreateFontCursor(dpy,XC_watch);
    XDefineCursor(dpy,win,cursor_watch);
    MyXSelectInput (dpy, win, (int) (ExposureMask
				       | KeyPressMask
				       | KeyReleaseMask
				       | ButtonPressMask
				       | ButtonReleaseMask
				       /*                               | ResizeRedirectMask   */
				       | StructureNotifyMask)
		  ); 
  font_info = XLoadQueryFont(dpy, "6x9");
  //XSetFont(dpy, gc, font_info->fid);
  {unsigned int first, last; 	   
  int id = font_info->fid;
    first = font_info->min_char_or_byte2;
    last = font_info->max_char_or_byte2;     
    fontList = glGenLists(last+1);
    if (fontList == 0) {
        printf ("out of display lists\n");
    exit (1);
    }
    glXUseXFont(id, first, last-first+1, fontList+first);     
    }
}    
    
#endif
#ifdef WGL
    a faire 
#endif
	carre = aspx == aspy;
	lacouleur = getcolor();
	nbcolor= 256; 
	ncolortable =0;
	LastColor=2;// En couleur pas defaul
	colortable=0;
	SetColorTable(2+6);

	INITGRAPH = 1;
    gluOrtho2D(0.0, height,0,width);
    glLineWidth(1);
	
    // cout <<" End Initgraphique\n";
}


void SetColorTable(int nb)
{
  if(!INITGRAPH) return;
   if (ncolortable == nb) return;// optim
   if (nbcolor && nb>2) 
     { 
       if(colortable) delete [] colortable;
       colortable = new RGBColor[nb];
       ncolortable = nb;
       if(LastColor>1) LastColor=nb-1;
       int k=0;
       colortable[k].red=65534;
       colortable[k].green=65534;
       colortable[k].blue=65534;
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
           //cout << "\t\t" << i0 << " " << j0 << " " << j1 << " " << k0 << " " << k1  << " "<<kk<<endl;
         	if(kk<=0) kk=1;
         // throwassert(kk);
         
          if (! grey)
           {
           colortable[k].red   = (cube6[j1][0]*k0+cube6[j0][0]*k1)/kk;
           colortable[k].green = (cube6[j1][1]*k0+cube6[j0][1]*k1)/kk;
           colortable[k].blue  = (cube6[j1][2]*k0+cube6[j0][2]*k1)/kk;
           }
          else 
           {
           kk=nb-1;
           k1 =  i0;
           k0 = nb - i0 -1;
           colortable[k].red   = (grey6[0][0]*k0+grey6[1][0]*k1)/kk;
           colortable[k].green = (grey6[0][1]*k0+grey6[1][1]*k1)/kk;
           colortable[k].blue  = (grey6[0][2]*k0+grey6[1][2]*k1)/kk;
           }
         
       /*    colortable[k].red   = (cube6[j1][0]*k0+cube6[j0][0]*k1)/kk;
           colortable[k].green = (cube6[j1][1]*k0+cube6[j0][1]*k1)/kk;
           colortable[k].blue  = (cube6[j1][2]*k0+cube6[j0][2]*k1)/kk;*/ 
           throwassert(k<ncolortable);
           
          }
  /*    for (k=0;k<ncolortable;k++)
           cout << " color"  << k 
                <<" r = " << colortable[k].red 
                <<" g = " << colortable[k].green
                <<" b = " << colortable[k].blue << endl;
  */    
         
       }
     else 
      ncolortable  =0;
}
void closegraphique(void)
{
  if(INITGRAPH) 
    {
    
 	DeleteFontGL (fontList);
#ifdef AGL 	
	aglSetCurrentContext (NULL);
	aglSetDrawable (ctx, NULL);
	aglDestroyContext (ctx);
    DisposeWindow(grafWindow0);
#endif
#ifdef XGL
      XUnloadFont(dpy, font_info->fid);
//      XFreeGC(dpy, gc);
      XCloseDisplay(dpy);
#endif
#ifdef WGL
  a faire 
#endif
      closePS();    
     delete [] colortable;colortable=0;
    }
  INITGRAPH=0;
}
void showgraphic()
{
#ifdef AGL
  if (grafWindow0 != FrontWindow())
   { 
	ShowWindow(grafWindow0);
	BringToFront(grafWindow0);
	SelectWindow(grafWindow0);
	SetPortWindowPort(grafWindow0); }
	GetPort(&grafPort0);
#endif
	
}

void reffecran(void)
{
  if(!INITGRAPH) return;
  		glClearColor(1.f, 1.f, 1.f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT); /* Clear buffer */

 
}

int getcolor(void)
{ return lacouleur;
}

void putpixel(int ix,int iy, int couleur)
{
//  if (ncolortable>3 && couleur < ncolortable && couleur >=0 ) 
//    SetCPixel(ix,iy,colortable+couleur);
//  DrawCStringGL ((char*) glGetString (GL_VENDOR), fontList);
  glBegin(GL_POINTS);
  glVertex2i(ix, height-iy);
  glEnd();

}

 void plotstring(const char *s)
{ 
// cout << "plotstring" << s << endl;
int lx=0,l = strlen(s);
 DrawCStringGL(s,fontList);
#ifdef XGL
  lx = XTextWidth( font_info,s,l);
#endif
 if(psfile) fprintf(psfile,"(%s) S\n",s);
 currx += lx;
} 

int LaCouleur(){return lacouleur;}

void couleur(int c)
{ 
  if ( lacouleur == c) // small optim
    return;
  c= c > LastColor ? 1 : c; // c=Min(c,LastColor); pour noir et blanc
  lacouleur =c;
    float r=1,g=1,b=1;
      if (c>=0 && c < ncolortable)
	{
	  r =  (float) colortable[c].red /65535.F;
	  g =  (float) colortable[c].green /65535.F;
	  b =  (float) colortable[c].blue /65535.F;
	 }
     else if (c!=0)
      r=g=b=0;
 	glColor4f (r,g,b,1.);
    if (psfile)
     fprintf(psfile,"%.3f %.3f %.3f C\n",r,g,b);
   
}

int InRecScreen(float x1, float y1,float x2, float y2)
{  
  float xi = Min(x1,x2),xa=Max(x1,x2);
  float yi = Min(y1,y2),ya=Max(y1,y2);
  return (xa >= rxmin) && (xi <= rxmax) && (ya >= rymin) && (yi <= rymax);
}
int InPtScreen( float x, float y)
{
  return (x >= rxmin) && (x <= rxmax) && (y >= rymin) && (y <= rymax);
}

void penthickness(int pepais)
{
//  PenSize(pepais,pepais);
  glLineWidth(pepais);
  if (psfile) fprintf(psfile,"%d setlinewidth\n",pepais);
}
void cadre(float xmin,float xmax,float ymin,float ymax)
{
  rxmin = xmin;
  rxmax = xmax;
  rymin = ymin;
  rymax = ymax;
  echx = aspx / (xmax - xmin);
  echy = aspy / (ymax - ymin);
}
void getcadre(float &xmin,float &xmax,float &ymin,float &ymax)
{
  xmin = rxmin;
  xmax = rxmax;
  ymin = rymin;
  ymax = rymax;

}

void cadreortho(float centrex, float centrey, float rayon)
{
  float xasp,yasp, getmaxx, getmaxy;
  
	getmaxx = xasp =aspx;	getmaxy = yasp = aspy;
	
  if (getmaxx * (float)xasp > getmaxy * (float)yasp)
  {
    rymin = centrey - rayon;
    rymax = centrey + rayon;
    echy= getmaxy / (2 * rayon);
    echx= (echy * xasp) / yasp;
    rxmin= centrex - getmaxx / (2 * echx);
    rxmax= centrex + getmaxx / (2 * echx);
  }
  else
  {
    rxmin = centrex - rayon;
    rxmax = centrex + rayon;
    echx = getmaxx / (2 * rayon);
    echy = (echx * yasp) / xasp;
    rymin = centrey - getmaxy / (2 * echy);
    rymax = centrey + getmaxy / (2 * echy);
  }
 // cout << "cadreortho\n";
}

int scalx(float x)
{
  return int((x - rxmin) * echx);
}

int scaly(float y)
{
  return int((rymax - y) * echy);
}

float scali(int i)
{
  return i/echx  + rxmin;
}
float scalj(int j)
{
  return -j/echy  + rymax;
}

void pointe(float x, float y)
{
  int newx = scalx(x), newy = scaly(y);
  putpixel(newx, newy, lacouleur);
  if (psfile) 
   fprintf(psfile,"%d %d P\n", newx, height-newy);
  
}

void rmoveto(float x, float y)
{
  int newx = scalx(x), newy = scaly(y);
 // MoveTo(newx,newy);
  if (psfile) 
   fprintf(psfile,"%d %d M\n", newx, height-newy);
  currx = newx; curry = newy;
  
}

void rlineto(float x, float y)
{
  int newx = scalx(x), newy = scaly(y);
  glBegin(GL_LINES);
  glVertex2i(currx, height-curry);
  glVertex2i(newx, height-newy);
  glEnd();
   if (psfile) 
    fprintf(psfile,"%d %d L\n", newx,height-newy);
  currx = newx; curry = newy;
  
}


void fillpoly(int n, float *poly)
{
  glBegin(GL_POLYGON);
  for (int i=0;i<n;i++)
    glVertex2i(scalx(poly[2*i]),height-scaly( poly[2*i+1]));
  glEnd();
    if (psfile) 
    {
     fprintf(psfile,"bF ");
     for (int i=0;i<n;i++)
      fprintf(psfile,"%d %d ", scalx(poly[2*i]),height-scaly( poly[2*i+1]));
     fprintf(psfile,"eF\n");
    }
}






int execute(const char* what)
{
  system(what);
  return 1; // error
}

#ifdef AGL
int DoMouseDown (int windowPart, WindowPtr whichWindow, EventRecord *myEvent)
{
  int wasactive;
	switch (windowPart) {
		case inGoAway:
			if (ours(whichWindow))
				if (TrackGoAway(whichWindow, myEvent->where))
					{ 
					   closegraphique();
					   cout << "Fin (fermeture fenetre graphique) " <<endl;
					   exit(0);
					   //HideWindow(whichWindow);  
					} 
			break;

		case inZoomIn:
			if (ours(whichWindow))
				{
					//SetCursor(&waitCursor);
					 SetPortWindowPort(whichWindow);
	                                 GetPort(&grafPort0);
					 reffecran();
					//EraseRect(&(whichWindow->portRect));
					ZoomWindow(whichWindow, inZoomIn, true);
					InitCursor();
				}
			break;

		case inZoomOut:
/*			if (ours(whichWindow))
				{
					SetCursor(&waitCursor); SetPort(whichWindow);
					EraseRect(&(whichWindow->portRect));
					ZoomWindow(whichWindow, inZoomOut, true);
					if(whichWindow == editWindow) 
						 MyZoomWindow(whichWindow);
					InitCursor();
				}*/
			break;

		case inMenuBar:
//			return(DoCommand(MenuSelect(myEvent->where)));
            break;

		case inSysWindow:
			//SystemClick(myEvent, whichWindow);
			break;

		case inDrag:
			if (ours(whichWindow))
				{
					SetPortWindowPort(whichWindow);
				        GetPort(&grafPort0);

		//			DragWindow(whichWindow, myEvent->where, &dragRect);
				}
			break;

		case inGrow:
			//if (ours(whichWindow))
			//	{MyGrowWindow(whichWindow, myEvent->where);}
			break;

		case inContent:
			wasactive = (whichWindow == FrontWindow()); 
	     if(!wasactive) { SelectWindow(whichWindow);
	//	    if (ours(whichWindow) && MacReDraw ) (* MacReDraw)();
		   }
		  else if (ours(whichWindow))
			{ SetPortWindowPort(whichWindow);	GetPort(&grafPort0);

			   while (Button()) ;
			   return 0;
			}
			break;
	}
return 1;
}
char HandleEvent(EventRecord	&	myEvent) 
{
  //  cout << "HandleEvent\n";
	WindowPtr		whichWindow=NULL;
	short			windowPart;
     
     char char1=0;
	   switch (myEvent.what) {
		case mouseDown:
  		windowPart = FindWindow(myEvent.where, &whichWindow);
	    if( DoMouseDown(windowPart, whichWindow, &myEvent) ==0) 
	      char1=  251;
	    break; 

//
//
//		case keyDown:
		case keyUp:
		case autoKey: 
			{

			 windowPart = FindWindow(myEvent.where, &whichWindow);
			if((whichWindow==grafWindow0) /* && (inContent == windowPart)*/)
				{ if  (grafWindow0 !=  FrontWindow()) { 
				     SelectWindow(whichWindow); 
		             SetPortWindowPort(whichWindow);	GetPort(&grafPort0);

		             }
				   char1 = (myEvent.message & 127L);
				  
				}
			break;}
		
	   case updateEvt:
	   /* 
	     if (ours((WindowPtr) myEvent.message)) {
	       	BeginUpdate((WindowPtr) myEvent.message);
		    EndUpdate((WindowPtr) myEvent.message);
	     } */
	    break;
   // cout << "End HandleEvent" << int(char1) << endl;

}
  return char1;
}


#endif

char Getijc(int & x,int & y)
{   
   char char1=0;

#ifdef AGL 
viderbuff();
   showgraphic();
    EventRecord		myEvent;
	int flag=1;
//	HLock( (Handle) WatchCurseur);
//	SetCursor(*CrossCurseur);
//	HUnlock( (Handle) WatchCurseur);
   SelectWindow(grafWindow0);
   while (char1==0) {
	if (GetNextEvent(everyEvent, &myEvent) /* ,OxFFFFFFFF,h)*/) 
	  char1=HandleEvent(myEvent);
	 }
   GlobalToLocal( & myEvent.where);
   x = myEvent.where.h;
   y = myEvent.where.v;
#endif
#ifdef XGL
  XEvent event;
  int flag,nb;
  XComposeStatus status;
  char buffer[20];
  KeySym keysym;   /*  incidence */
  XDefineCursor(dpy,win,cursor_arrow);
  flag=0;
  while (!flag)
  { XNextEvent(dpy, &event);
    if(event.type == ButtonRelease) 
    { x = event.xbutton.x;
      y = event.xbutton.y; 
      if      (event.xbutton.button == Button1) char1=shift?248:251;
      else if (event.xbutton.button == Button2) char1=shift?249:252;
      else                                      char1=shift?250:253; 
      //     printf(" mouse release %d\n",(int) char1);
      flag=1;
    }
    else if(event.type == KeyPress)
    { x = event.xkey.x;
      y = event.xkey.y; 
      char1= event.xkey.keycode ;
       keysym=0;
       nb=XLookupString(&event.xkey,buffer,20,&keysym,&status);

/*        printf("nb= %d keysym= %d buffer=",nb,keysym);
/*        for(i=0;i<20;i++)
/*         printf(" %d ",(int)buffer[i]);
/*        printf("\n");
*/

/*       voir    /usr/include/X11/keysymdef.h + ap_keysym */

       if (nb != 0) 
         {char1 = buffer[0];
          flag= 1; 
         }
       else
         {
/*          if     (IsFunctionKey(keysym))     printf("function down\n");
          else if(IsModifierKey(keysym))     printf("modifier down\n");
          else if(IsKeypadKey(keysym))       printf(" keypad down\n");
          else if(IsMiscFunctionKey(keysym)) printf(" misc function down\n");
          else if(IsPFKey(keysym))           printf(" PF key down\n");
*/
#ifdef XK_MISCELLANY	      

          switch(keysym) 
            {
/* Cursor control & motion */
	      /*
            case XK_Left :
              flag = 1;
              char1 = call(keyboa).curs_left;
              break;
            case XK_Up :
              flag = 1;
              char1 = call(keyboa).curs_up;
              break;
            case XK_Right :
              flag = 1;
              char1 = call(keyboa).curs_right;
              break;
            case XK_Down :
              flag = 1;
              char1 = call(keyboa).curs_down;
              break;
            case XK_Next :
              flag = 1;
              char1 = call(keyboa).pad_down;
              break;
            case XK_Prior :
              flag = 1;
              char1 = call(keyboa).pad_up;
              break;
            case XK_End :
              flag = 1;
              char1 = call(keyboa).marg_right;
              break;
            case XK_Begin :
              flag = 1;
              char1 = call(keyboa).marg_left;
              break;
	      */
/* Misc Functions */ 
	      /* 

            case XK_Select :
              flag = 1;
              char1 = call(keyboa).mark;
              break; */
/*
            case XK_Print :
              flag = 1;
              char1 = ;
              break;  
            case XK_Execute :
              flag = 1;
              char1 = ;
              break;  
            case XK_Insert :
              flag = 1;
              char1 = ;
              break;

            case XK_Undo :
              flag = 1;
              char1 = call(keyboa).undo;
              break;

            case XK_Redo :
              flag = 1;
              char1 = ;
              break;
            case XK_Menu :
              flag = 1;
              char1 = ;
              break;
            case XK_Find :
              flag = 1;
              char1 = ;
              break;


            case XK_Cancel :
              flag = 1;
              char1 = call(keyboa).line_del;
              break;
            case XK_Help :
              flag = 1;
              char1 = call(keyboa).help;
              break;

            case XK_Break :
              flag = 1;
              char1 = ;
              break;
            case XK_Mode_switch :
              flag = 1;
              char1 = ;
              break;
            case XK_script_switch :
              flag = 1;
              char1 = ;
              break;
            case XK_Num_Lock :
              flag = 1;
              char1 = ;
              break;

            case XK_F1 :
              flag = 1;
              char1 = shift ? call(keyboa).sfunct1 : call(keyboa).funct1 ;
              break;
            case XK_F2 :
              flag = 1;
              char1 = shift ? call(keyboa).sfunct2 : call(keyboa).funct2 ;
              break;
            case XK_F3 :
              flag = 1;
              char1 = shift ? call(keyboa).sfunct3 : call(keyboa).funct3 ;
              break;
            case XK_F4 :
              flag = 1;
              char1 = shift ? call(keyboa).sfunct4 : call(keyboa).funct4 ;
              break;
            case XK_F5 :
              flag = 1;
              char1 = shift ? call(keyboa).sfunct5 : call(keyboa).funct5 ;
              break;
            case XK_F6 :
              flag = 1;
              char1 = shift ? call(keyboa).sfunct6 : call(keyboa).funct6 ;
              break;
            case XK_F7 :
              flag = 1;
              char1 = shift ? call(keyboa).sfunct7 : call(keyboa).funct7 ;
              break;
            case XK_F8 :
              flag = 1;
              char1 = shift ? call(keyboa).sfunct8 : call(keyboa).funct8 ;
              break;
            case XK_F9 :
              flag = 1;
              char1 = shift ? call(keyboa).sfunct9 : call(keyboa).funct9 ;
              break;
            case XK_F10 :
              flag = 1;
              char1 = shift ? call(keyboa).sfunct10 : call(keyboa).funct10 ;
              break;
            case XK_F11 :
              flag = 1;
              char1 = shift ? call(keyboa).sfunct11 : call(keyboa).funct11 ;
              break;
            case XK_F12 :
              flag = 1;
              char1 = shift ? call(keyboa).sfunct12 : call(keyboa).funct12 ;
              break;
	      */
            case XK_Shift_L :
              shift=1;
              break;
            case XK_Shift_R :
              shift=1;
              break;
           case XK_Control_L :
              control=1;
              break;
            case XK_Control_R :
              control=1;
              break;
            case XK_Caps_Lock :
              shiftlock = 1 ;
              break;
            case XK_Shift_Lock :
              shiftlock = 1 ;
              break;
            case XK_Meta_L :
              alt=1;
              break;
            case X<K_Meta_R :
              alt=1;
              break;
            case XK_Alt_L :
              alt=1;
              break;
            case XK_Alt_R :
              alt=1;
              break;
            } /* end switch */
#endif              
         }
    }
    else if(event.type == KeyRelease)
    { x = event.xkey.x;
      y = event.xkey.y; 
      char1= event.xkey.keycode ;
       keysym=0;
       nb=XLookupString(&event.xkey,buffer,20,&keysym,&status);
/*          if     (IsFunctionKey(keysym))     printf("function up\n");
          else if(IsModifierKey(keysym))     printf("modifier up\n");
          else if(IsKeypadKey(keysym))       printf(" keypad up\n");
          else if(IsMiscFunctionKey(keysym)) printf(" misc function up\n");
          else if(IsPFKey(keysym))           printf(" PF key up\n");
*/
       if (nb == 0) 
         {
#ifdef XK_MISCELLANY	      
          switch(keysym)
            {
            
            case XK_Shift_L :
              shift=0;
              break;
            case XK_Shift_R :
              shift=0;
              break;
           case XK_Control_L :
              control=0;
              break;
            case XK_Control_R :
              control=0;
              break;
            case XK_Caps_Lock :
              shiftlock = 0 ;
              break;
            case XK_Shift_Lock :
              shiftlock = 0 ;
              break;
            case XK_Meta_L :
              alt=0;
              break;
            case XK_Meta_R :
              alt=0;
              break;
            case XK_Alt_L :
              alt=0;
              break;
            case XK_Alt_R :
              alt=0;
              break;
              
            } /* end switch */
#endif	      

         }
    }
  }
  XDefineCursor(dpy,win,cursor_watch);
  XFlush(dpy);
#endif
#ifdef WGL
  a faire 
#endif
  return char1;
    
    

}


char Getxyc(float &x,float &y)
{ 
  char c;
  int i,j;
 // cout << "getxyc \n";
  c = Getijc( i,j);
  x = scali(i);
  y = scalj(j);
//  cout << "getxyc out \n";
  return c;
}

void  viderbuff(){
     glFinish();
#ifdef AGL     
     aglSwapBuffers (ctx); // send swap command
#endif
#ifdef XGL     
     glXSwapBuffers (dpy,win); // send swap command
#endif

}


void rattente(int waitm)
{ int i,j;
 char   c=0;
 if(waitm)  c = Getijc( i,j);
 if ( c == 3) {cout << "rattente: ^c => abort " << endl;closegraphique();exit(1);}// ^c  => exit
/*    you may prefer to use carriage return to move to the next graph */
/*	 getc(stdin);
*/
// if(waitm) while(!Button()){ };
}

void GetSizeScreen(int & ix,int &iy);
void GetSizeScreen(int & ix,int &iy)
{
  	ix = width ;
  	iy = height;
}



void openPS(const char *filename )
{ 
  char ffff[32];
  int count=0;
  if(psfile_save) closePS();
  time_t t_loc;
  float s=0.5;
  const int shiftx=50,shifty=50;
 // char  username[10];
  time(&t_loc);
  bool notfound;
  if( !filename) 
   do {
      struct stat buf;
      sprintf(ffff,"rgraph_%.3d.ps",count++);
      volatile int r= stat(ffff,&buf) ;
      notfound = r !=0;
      if(count>1000) break;
    } while ( !notfound );
   
  psfile=fopen(filename?filename:ffff,"w");
   
  if(psfile==0) {printf("Erreur %s errno %d\n",filename?filename:ffff,errno);exit(1);}
  if(psfile) {
  fprintf(psfile,"%%!PS-Adobe-2.0 EPSF-2.0\n%%%%Creator: %s\n%%%%Title: FreeFem++\n","user");
  fprintf(psfile,"%%%%CreationDate: %s",ctime(&t_loc));
  fprintf(psfile,"%%%%Pages: 1\n");
  fprintf(psfile,"%%%%BoundingBox:       %d %d %d %d\n",shiftx,shifty,int(shiftx+width*s),int(shifty+height*s));
  fprintf(psfile,"%%%%EndComments\n");
  fprintf(psfile," /L {  lineto currentpoint stroke newpath moveto} def\n");
  fprintf(psfile," /M {  moveto } def\n");
  fprintf(psfile," /C {setrgbcolor} def\n");
  fprintf(psfile," /rec {newpath 4 copy 8 1 roll moveto 3 -1 roll lineto 4 2 roll exch lineto lineto closepath} def\n");
  fprintf(psfile," %d %d  translate \n",shiftx,shifty);
  fprintf(psfile," %f %f  scale \n",s,s);
  fprintf(psfile," 0 %d 0 %d rec clip newpath\n",int(width),int(height));
  fprintf(psfile," /Helvetica findfont 10 scalefont setfont\n");
  fprintf(psfile," /S { show} def\n");
  fprintf(psfile," /bF  { mark} def \n");
  fprintf(psfile," /eF {newpath moveto counttomark 2 idiv {lineto} repeat closepath fill cleartomark} def\n");
  fprintf(psfile," /P { /yy exch def /xx exch def   xx xx 1 add yy yy 1 add  rec  fill } def\n");

  fprintf(psfile," 1 setlinewidth\n");
  psfile_save=psfile;
  }
}
void closePS(void)
{
  if(psfile_save)   {
    fprintf(psfile_save,"showpage\n");
    fclose(psfile_save);
    }
    
  psfile=0;
  psfile_save=0;
  
}

  void Commentaire(const char * c)  
  {
  if(psfile)   {
    fprintf(psfile,"%% %s\n",c);
   }
  };
  void NoirEtBlanc(int NB)
  {
    if(NB) LastColor=1;
    else LastColor=ncolortable?ncolortable:2;
  }
 
  void MettreDansPostScript(int in)
   {
     if(in)  psfile=psfile_save;     
     else    psfile=0;
   }

static void     FillRect(float x0,float y0, float x1, float y1)
 {
     float r[8];
     r[0]=x0;r[1]=y0;
     r[2]=x1;r[3]=y0;
     r[4]=x1;r[5]=y1;
     r[6]=x0;r[7]=y1;
     fillpoly(4,r);
 }

float  GetHeigthFont()
{ 
 // FontInfo 	MyFontInfo;
 // GetFontInfo(&MyFontInfo);  
#ifdef XGL   				
  int dir,asc,desc,k;
  XCharStruct overall;

  XTextExtents(font_info,"gML",3,&dir,&asc,&desc,&overall); 
  return (asc+desc)*(0.9/echy);
#else  
 int interligne = 9;// MyFontInfo.ascent + MyFontInfo.descent + MyFontInfo.leading;
 return interligne/echy;
#endif

}



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
  sprintf(buf,"%g",xf);
  couleur(1);
  plotstring(buf);

   return lineno;
}
 void ShowHelp(const char * s,int k)
{
  if(k) {
    MettreDansPostScript(0);
    couleur(1);
    float xmin,xmax,ymin,ymax;
    getcadre(xmin,xmax,ymin,ymax);
    rmoveto(xmin+(xmax-xmin)/100,ymax-(k)*(ymax-ymin)/30);
    plotstring(s);
    MettreDansPostScript(1);
       //  couleur(1);	
  }
}

class Grid;
  void setgrey(bool gg ){grey=gg;}
  int getgrey(){ return grey;}

void SaveMesh(Grid &t){}
void SavePlot(int D, Grid& t, double *f){}
void SavePlot(int D, Grid& t, float *f){}

