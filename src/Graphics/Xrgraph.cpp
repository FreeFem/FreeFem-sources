// ********** DO NOT REMOVE THIS BANNER **********
//
// SUMMARY: Bamg: Bidimensional Anisotrope Mesh Generator
// RELEASE: 0 
// USAGE  : You may copy freely these files and use it for    
//          teaching or research. These or part of these may   
//          not be sold or used for a commercial purpose with- 
//          out our consent : fax (33) 1 39 63 55 14       
//
// AUTHOR:   D. Bernardi, F. Hecht,  O. Pironneau ,    Y. Darmaillac                      
// ORG    :  INRIA
// E-MAIL :   Frederic.Hecht@Inria.fr   
//
// ORIG-DATE:     Dec 97

#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <cassert>
using namespace std;

#define MAXSHORT 0xFFFF

#ifdef HPPA
#ifndef __GNUC__
typedef char *caddr_t;
#endif
#endif
#ifdef __MWERKS__
#include <Xlib.h>
#include <Xutil.h>
#include <Xos.h>
#include <Xatom.h>
#include <keysym.h>
#include <cursorfont.h>
#else
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <X11/Xatom.h>
#include <X11/keysym.h>
#include <X11/cursorfont.h>
#endif
#include "rgraph.hpp"

#ifdef macintoshxx
#include <ConditionalMacros.h>
#include <unix.h>
#else
#include <sys/stat.h>
#endif

template<class T> inline T Min (const T &a,const T &b){return a < b ? a : b;}
template<class T> inline T Max (const T &a,const T & b){return a > b ? a : b;}


static  long  cube6[7][3] ={ { 65535,32000,32000},{ 65535, 65535,0},{0, 65535,0},{0, 65535, 65535},{0,0, 65535}
     , { 65535,0, 65535},{ 32000,0,0} }; 
static  long grey6[2][3] ={ {65534,65534,65534},{0,0,0} }; 

static FILE *psfile = 0;
static FILE *psfile_save = 0;
static bool grey=false;
static int LastColor=2;  //  pour est en couleur par defaut

const float fMinPixel = -32000;
const float fMaxPixel = +32000;
#define reel float 
static  Display *display;
static  Window win;
static  XSizeHints size_hints;
// static  XEvent report;
static int ncolortable,fcolor;
static XColor *colortable;
static  GC gc;
static  XFontStruct *font_info;
static int shift, control,shiftlock,alt;
static reel echx,echy,rxmin,rxmax,rymin,rymax;
static int  lacouleur,screen, width, height, currx, curry;
static unsigned long  background,foreground;
static Cursor cursor_watch,cursor_arrow;
static long NbErrX11 =0;
Colormap color_map,color_map_sys;
#define call(i) i
static  Visual *visual;
static int INITGRAPH=0;
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
extern int mymain(int argc,char **argv);
int main (int argc, char **argv)
{
 return mymain(argc,argv);
}

#endif

void message(char *s)
{  printf("%s	\n",s);}

void erreur(char *s)
{ message(s); exit(0);}

void *safecalloc(size_t nb, size_t  size)
{
  void* p=NULL;
  p = calloc(nb, size);
  if (p == NULL) printf("Run out of Memory!\n");
  return p;
}

void safefree(void** f)
{
  if(*f){ free((char*) *f); *f=NULL;}
}

void rflush()
{
  XEvent report;
  XNextEvent(display, &report);
  if (report.type == Expose)
    while (XCheckTypedEvent(display, Expose, &report));
  
  XFlush(display);
}

int xerror  (Display *display,XErrorEvent * myerr)
{
  if (NbErrX11++<10) {
  char msg[80];
  XGetErrorText(display, myerr->error_code, msg, 80);
  fprintf(stderr, "Error code %s\n", msg);}
  return 0;
}

/*
void xerror()
{
  fprintf(stderr, "Probleme avec X-Windows\n");
  assert(0);
}
*/
void xerrorio()
{
  
  fprintf(stderr, "Fatal erreur avec X-Windows\n");
  assert(0);
  exit(2);
}
void  MyXSelectInput(Display * dpy,Window w,int  mask)

{
  XSetWindowAttributes attributes;
  attributes.event_mask = mask;
  XChangeWindowAttributes(dpy, w, CWEventMask, &attributes);
}
int LaCouleur() {return  lacouleur;}

void couleur(int c)
{ 
  if ( lacouleur == c) // small optim
    return;

  c= c > LastColor ? 1 : c; // c=Min(c,LastColor); pour noir et blanc
 lacouleur = c;
 if (colortable)
   { 
     if (c>=0 && c < ncolortable)
       XSetForeground(display,gc,colortable[c].pixel);
     else 
       XSetForeground(display,gc,foreground);
   }
else
  if ( c == 0 )
    XSetForeground(display,gc,background);
  else
    XSetForeground(display,gc,foreground);
 if (psfile)
  {
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
    
    fprintf(psfile,"%.3f %.3f %.3f C\n",r,g,b);
  }
}

void SetColorTable(int nb)
{
  int i;
   if (fcolor  && nb>2 && nb < 256) 
     { 
       nb = Max(nb,8);
       if (ncolortable == nb) return;// optim
       if(colortable) delete [] colortable;
       colortable = new XColor[nb];
       ncolortable = nb;
       if(LastColor>1) LastColor=nb-1;
       int k=0;
       colortable[k].pixel=k;
       colortable[k].red= 65535;
       colortable[k].green= 65535;
       colortable[k].blue= 65535;
       colortable[k].flags = DoRed | DoGreen | DoBlue;
       background=k;
       k++;
       colortable[k].pixel=k;
       colortable[k].red=0;
       colortable[k].green=0;
       colortable[k].blue=0;
       colortable[k].flags = DoRed | DoGreen | DoBlue;
       foreground=k;
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
	   colortable[k].pixel=  (unsigned long) k;
	   colortable[k].flags = DoRed | DoGreen | DoBlue;
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
       if (visual->c_class != TrueColor)
	 {
	   // cout << "XStoreColors( not TrueColor)" << ncolortable << " "  <<
	     XStoreColors (display, color_map, colortable, ncolortable) ;
	       //	<< endl;
	 }
       else 
	 {
	   // cout << "XAllocColor (TrueColor)" << endl; 
	   for (i=0;i<ncolortable;i++)
	     XAllocColor(display, color_map, colortable+i );
	   background =colortable[background].pixel;
	   foreground =colortable[foreground].pixel;
	 }
       if (win) {
	 XGCValues gcvalues;
	 gcvalues.foreground = foreground;
	 gcvalues.background = background;
	 gc = XCreateGC(display, win, GCForeground | GCBackground , &gcvalues);
       }
     }


// a faire 
}
void FlushEvent()
{
  XEvent event;
  while (XPending(display))
   XNextEvent(display, &event);
} 

void initgraphique()
{
 int ddd;
  atexit(myend);
  win=0;
  XSetWindowAttributes attributes;
  NbErrX11=0;
  XGCValues gcvalues;
  XEvent report;
  display = XOpenDisplay(NULL);
  if (!display) 
    {
      cerr << " Erreur openning  display " << endl;
      exit(2);
    }
  colortable=0;
  ncolortable=0;
  LastColor=2;// En couleur pas defaul
  //  modif FH
  Display *dpy=display;
  // Colormap color_map,color_map_sys;
  visual  = DefaultVisual(display, DefaultScreen(display));
  fcolor=0;  /* pas couleur */
  int fstereo=0; /* non */
  int nbplans=visual->bits_per_rgb;
  color_map_sys = DefaultColormap (display, DefaultScreen (display)); 
  color_map = color_map_sys;
  foreground= BlackPixel(display, screen);
  background= WhitePixel(display, screen);

  switch (visual->c_class)
	 {
	 case GrayScale:  {break;}
	 case PseudoColor: 
	   { 
	     cout << " PseudoColor  nbcolor =" << visual->map_entries << endl;
	     color_map= XCreateColormap (display, RootWindow (display, DefaultScreen (display)),
					 visual,AllocAll);
	     // copy the def color map 
	     for (int i=0;i<visual->map_entries;i++)
	       {  XColor colorcell_defs;
	          colorcell_defs.pixel = (unsigned long) i;
		  XQueryColor (display, color_map_sys, &colorcell_defs);
		  XStoreColor (display, color_map, &colorcell_defs);
	       }
	     fcolor=1;
	     
	     SetColorTable(8); // set 
	     break;
	   }
	 case DirectColor:
	   {
	     cout << " DirectColor " << endl;
	     fcolor=1;
	     SetColorTable(8); // set 
	     break; 
	   }
	 case TrueColor  : 
	   {
	     cout << " TrueColor " << endl;
	     fcolor=1;	     
	     SetColorTable(8); // set 
	     break;
	   }
	 } 
  font_info = XLoadQueryFont(display, "6x9");
  XSetErrorHandler((XErrorHandler)xerror);
  XSetIOErrorHandler((XIOErrorHandler)xerrorio);
  screen = DefaultScreen(display);
  width = DisplayWidth(display, screen);
  height = DisplayHeight(display, screen);
  ddd = width < height ?  width : height;
  width = ddd*8/10;
  height =  ddd*8/10;


  attributes.background_pixel = background;
  attributes.border_pixel     = foreground;
  attributes.backing_store    = Always;
  attributes.colormap         = color_map;

  win = XCreateWindow(display, RootWindow(display, DefaultScreen(display)),
		      50, 80, width, height,4,
		      CopyFromParent, InputOutput, visual,
		      CWBackPixel | CWBorderPixel | CWBackingStore | CWColormap,
		      &attributes);
  const char * title("FreeFem 1.0");
  XChangeProperty(display, win, XA_WM_NAME, XA_STRING, 8
                  , PropModeReplace,(const unsigned char *)  title , strlen(title)); 


  
         
  gcvalues.foreground = foreground;
  gcvalues.background = background;
  gcvalues.function   = GXcopy    ;
  gc = XCreateGC(display, win, GCForeground | GCBackground | GCFunction, &gcvalues);
  
  XSetFillRule(display,gc,WindingRule);
  
  
  // win = XCreateSimpleWindow(display, RootWindow(display, screen), 50, 80, width, height, 4,
  // foreground,background);
  cursor_arrow = XCreateFontCursor(display,XC_arrow);
  cursor_watch = XCreateFontCursor(display,XC_watch);
    
  size_hints.flags = PPosition | PSize;
  size_hints.x = 0;
  size_hints.y = 0;
  size_hints.width = width;
  size_hints.height = height;

  XSetFont(display, gc, font_info->fid);
  XSetForeground(display, gc, foreground); 
  XMapWindow(display, win);
  MyXSelectInput (display, win, (int) (ExposureMask
				       | KeyPressMask
				       | KeyReleaseMask
				       | ButtonPressMask
				       | ButtonReleaseMask
				       /*                               | ResizeRedirectMask   */
				       | StructureNotifyMask)
		  ); 
  
  // do XNextEvent(display, &report); while (report.type != Expose);
  XDefineCursor(display,win,cursor_watch);
  XFlush(display); 
  INITGRAPH = 1;

}

void closegraphique()
{
  if (INITGRAPH)
    {
      INITGRAPH = 0;
      XUnloadFont(display, font_info->fid);
      XFreeGC(display, gc);
      XCloseDisplay(display);
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

void pointe(reel x, reel y)
{
  XDrawPoint(display, win, gc, scalx(x), scaly(y));
}

void rmoveto(reel x, reel y)
{
  currx = scalx(x);
  curry = scaly(y);
}

void rlineto(reel x, reel y)
{
  int newx = scalx(x), newy = scaly(y);
  XDrawLine(display, win, gc, currx, curry, newx, newy);
  if (psfile)
    fprintf(psfile,"%d %d %d %d L\n",currx, height-curry, newx, height-newy);
  currx = newx; curry = newy;
/*   XFlush(display); */
}

void cadreortho(reel centrex, reel centrey, reel rayon)
{
  //  int xasp,yasp;

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
{ int lx,l = strlen(string);
  XDrawString(display, win, gc, currx, curry , string, l);
 lx = XTextWidth( font_info,string,l);
 if(psfile) fprintf(psfile,"(%s) %d %d  S\n",string,currx,height-curry);
 currx += lx;
}

void showgraphic()
{
}

void x11draw3(int * ptype)
{
  XGCValues gcvalues;
  int type;

  type=  *ptype;
  switch (type)
  {
    case 0  : {gcvalues.line_style = LineSolid;     break;}
    case 1  : {gcvalues.line_style = LineOnOffDash; break;}
    default : {gcvalues.line_style = LineDoubleDash;break;}
  }
  XChangeGC(display, gc, GCLineStyle, &gcvalues);

  if (psfile) 
    switch (type) {
    case 0  : {fprintf(psfile,"[] setdash\n");break;}
    case 1  : {fprintf(psfile,"[3]  setdash\n");break;}
    default : {fprintf(psfile,"[4 1] setdash\n");break;}
    }

}  

void penthickness(int pepais)
{
  XGCValues gcvalues;
  gcvalues.line_width = pepais;
  XChangeGC(display, gc, GCLineWidth, &gcvalues);
  if (psfile) fprintf(psfile,"%d setlinewidth\n",pepais);
}


void x11linsrn(int * x1,int * x2,int * y1,int * y2)
  //int *x1,*x2,*y1,*y2;
{   
  XDrawLine(display, win, gc, *x1, *x2, *y1, *y2); 
/*  call(viderbuff)(); */
}

   
void viderbuff()
{
  XRaiseWindow (display,win);
  XFlush(display); 
}



void cercle(reel centrex, reel centrey, reel rayon)
{
  int r = (int) (rayon * echx);
  XDrawArc(display, win, gc,
	   scalx(centrex) - r, scaly(centrey) - r, width, height, 0, 360 * 64);
  XFlush(display);
}
void reffecran()
{
 XClearWindow(display,win);
}

void fillpoly(int n, float *poly)
{
  int i;
  XPoint *poly0,polyloc[10];
  if(n<10)
    poly0=polyloc;
  else
    if(poly0= (XPoint *) malloc(n*sizeof(XPoint)), !poly)
  {
    fprintf(stderr, "Erreur d'allocation dans raffpoly\n");
    return;
  }
  for(i=0; i<n; i++)
  {
    poly0[i].x =scalx(poly[2*i]);
    poly0[i].y =scaly(poly[2*i+1]);
  }
  
  XFillPolygon(display, win, gc, poly0,n, Complex, CoordModeOrigin);

  if( poly0!=polyloc) free((char*)poly0);
   if (psfile) 
    {
     fprintf(psfile,"bF ");
     for (i=0;i<n;i++)
      fprintf(psfile,"%d %d ", scalx(poly[2*i]),height-scaly( poly[2*i+1]));
     fprintf(psfile,"eF\n"); 
    }
}

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


int execute (const char * str)
{ 
  return system(str);
}

char Getijc(int *x1,int *y1)
{ char char1;
  XEvent event;
  int flag,nb;
  XComposeStatus status;
  char buffer[20];
  KeySym keysym;   /*  incidence */
  XDefineCursor(display,win,cursor_arrow);
  flag=0;
  while (!flag)
  { XNextEvent(display, &event);
    if(event.type == ButtonRelease) 
    { *x1 = event.xbutton.x;
      *y1 = event.xbutton.y; 
      if      (event.xbutton.button == Button1) char1=shift?248:251;
      else if (event.xbutton.button == Button2) char1=shift?249:252;
      else                                      char1=shift?250:253; 
      //     printf(" mouse release %d\n",(int) char1);
      flag=1;
    }
    else if(event.type == KeyPress)
    { *x1 = event.xkey.x;
      *y1 = event.xkey.y; 
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
            case XK_Meta_R :
              alt=1;
              break;
            case XK_Alt_L :
              alt=1;
              break;
            case XK_Alt_R :
              alt=1;
              break;
            } /* end switch */
         }
    }
    else if(event.type == KeyRelease)
    { *x1 = event.xkey.x;
      *y1 = event.xkey.y; 
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

         }
    }
  }
  XDefineCursor(display,win,cursor_watch);
  XFlush(display);
  return char1;
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
void rattente(int waitm)
{ int i,j;
  XFlush(display);
  if (waitm)
      Getijc(&i,&j);
}
 void GetScreenSize(int &ix,int &iy)
{
  ix = width;
  iy = height;
}
void openPS(const char *filename )
{ 
  char ffff[32];
  int count=0;
  if(psfile_save) closePS();
  time_t t_loc;
  float s=0.5;
  char  username[10];
  /*if (!cuserid(username)) */ strcpy(username,"inconnue");
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
   

  const char *fps (filename?filename:ffff);

  
  psfile=fopen(fps,"w");
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
    fprintf(psfile," /Helvetica findfont 10 scalefont setfont\n");
    fprintf(psfile," /S {moveto show} def\n");
    fprintf(psfile," /bF  { mark} def \n");
    fprintf(psfile," /eF {newpath moveto counttomark 2 idiv {lineto} repeat closepath  fill cleartomark} def\n");
    fprintf(psfile," /P { /yy exch def /xx exch def   xx xx 1 add yy yy 1 add  rec  fill } def\n");
    fprintf(psfile," 1 setlinewidth\n");
  }
  else 
    cerr << " Err openning postscript file " << fps << endl;
}
void closePS(void)
{
  if(psfile_save) {
    fprintf(psfile_save,"showpage\n");
    fclose(psfile_save);
  }
  psfile_save=0;
  psfile=0;
}

 void coutmode(short i)  {}
// bof bof --- 
 float  GetHeigthFont()
{ 
  int dir,asc,desc,k;
  XCharStruct overall;

  XTextExtents(font_info,"gML",3,&dir,&asc,&desc,&overall); 
  return (asc+desc)*(0.9/echy);
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
     else   psfile=0;
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

  void setgrey(bool gg ){grey=gg;}
  int getgrey(){ return grey;}

class Grid;

void SaveMesh(Grid &t){}
void SavePlot(int D, Grid& t, double *f){}
void SavePlot(int D, Grid& t, float *f){}

