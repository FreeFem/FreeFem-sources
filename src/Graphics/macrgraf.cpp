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
//#define TARGET_API_MAC_CARBON 1
#include <MSLCarbonPrefix.h>
#include <sioux.h>
#include <SIOUXGlobals.h> //OP my hack
#include "error.hpp"
#include <fstream>
#include <iostream>
#include <stdlib.h>
using namespace std;

#include <errno.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#define fill thequikdrawfill
#include "rgraph.hpp"
#include <Windows.h>
#include <Fonts.h>
#include <SegLoad.h>
#include <Quickdraw.h>
//#include <StandardFile.h>
#include <Navigation.h>
#include <time.h>
#include <setjmp.h>
#include <time.h>
#include <unix.h>
#include <SIOUX.h> 
#undef fill

static FILE *psfile = 0;
static FILE *psfile_save = 0;

static bool grey=false;
int pStrCopy (StringPtr p1, char * p2);

#ifdef FREEFEM
//  pour imprimer la version   FH 
#define STRING(i) #i
#include <new.h>

jmp_buf environ;
static int  myenviron = 0;

TEHandle TESioux;

void out_of_memory ();
void NEW_HANDLER (void);
void compile(char *fname);
float scali(int i);
float scalj(int j);
void execute(char* what);
int DoMouseDown (int windowPart, WindowPtr whichWindow, EventRecord *myEvent);
char Getijc(int & x,int & y);
 

void out_of_memory ()
{
  cout << "FreeFEM error: operator new failed; not enough memory" << endl;
  if (myenviron)
   longjmp(environ,1);
  exit(2);
}

void NEW_HANDLER (void){  set_new_handler (&out_of_memory);}
#endif

#define	ours(w)		(w==grafWindow0)

template<class T> inline T Min (const T &a,const T &b){return a < b ? a : b;}
template<class T> inline T Max (const T &a,const T & b){return a > b ? a : b;}

static  int cube6[7][3] ={ {65534,0,0},{65534,65534,0},{0,65534,0},{0,65534,65534},{0,0,65534}
     , {65534,0,65534},{65534,0,0} }; 
static  int grey6[2][3] ={ {65534,65534,65534},{0,0,0} }; 

char errbuf[255];
static int INITGRAPH=0;
static float aspx, aspy, echx,echy,ech,rxmin,rxmax,rymin,rymax;
static int carre, lacouleur;
// static CWindowRecord wgRecord0;
static WindowPtr	 grafWindow0;
static GrafPtr          grafPort0;
static	Rect		boundsRect;
static int nbcolor;
static CursHandle  CrossCurseur ;
static CursHandle  WatchCurseur ;
static int ncolortable;
static int LastColor; // LastColor=1 => Noir et Blanc 
static  Pattern  white,black;

static int width,height;
static RGBColor * colortable;
int getcolor();
void putpixel(int ix,int iy, int couleur);
int scalx(float x);
int scaly(float y);
void thisexit();
void InitMac();
// --------------------------------------------------------------------------------------------------------------

// --------------------------------------------------------------------------------------------------------------
// APPLE EVENT SUPPORT ROUTINES
// --------------------------------------------------------------------------------------------------------------
static OSStatus	MissingParameterCheck(
	const AppleEvent 	*inputEvent)
/*
	This routine checks an input AppleEvent for the missing keyword.
	If the missing keyword is found, that means that some required
	parameters were missing (ie, an error). 
	
	However, if the missing keyword isn't found, that means that we aren't missing 
	any required parameters (that is to say, all REQUIRED parameters were supplied
	by the person who created the event).
	
	SOME DAY, THE ABOVE COMMENT WILL MAKE SENSE TO YOU.  IT STILL DOESN'T
	TO ME AND I WAS THE ONE WHO WROTE IT.
*/
{
	OSStatus		anErr;
	AEKeyword	missingKeyword;
	DescType	ignoredActualType;
	Size		ignoredActualSize;
	
	anErr = AEGetAttributePtr(
		inputEvent, 
		keyMissedKeywordAttr,
		typeWildCard,
		&ignoredActualType,
		(Ptr) &missingKeyword,
		sizeof(AEKeyword),
		&ignoredActualSize);
			
	if (anErr == noErr)
		anErr = errAEParamMissed;
	else
		if (anErr == errAEDescNotFound)
			anErr = noErr;
		
	return anErr;
	
} // MissingParameterCheck


static pascal OSErr	DoOpenApp(
	const AppleEvent 	*inputEvent,
	AppleEvent 	*outputEvent,
	SInt32		handlerRefCon)
{
/*
 #pragma unused (outputEvent, handlerRefCon)

	DoCommand(nil, cNew, 0, 0);
	
	// so that the initial document opens more quickly, we don't start
	// the threads until we get an OpenApp or OpenDocument AppleEvent
	if (gStarterThread != kNoThreadID)
		SetThreadState(gStarterThread, kReadyThreadState, gStarterThread);
*/	
	return(MissingParameterCheck(inputEvent));
	
} // DoOpenApp

// --------------------------------------------------------------------------------------------------------------
static pascal OSErr	DoReopenApp(
	const AppleEvent 	*inputEvent,
	AppleEvent 	*outputEvent,
	SInt32		handlerRefCon)
{
/*
#pragma unused (outputEvent, handlerRefCon)

	if (FrontWindow() == nil)
		DoCommand(nil, cNew, 0, 0);
*/	
	return(MissingParameterCheck(inputEvent));
	
} // DoReopenApp

// --------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------
static pascal OSStatus	DoOpenFile(
	const AppleEvent 	*inputEvent
	)	// nil == 0, zero length == print to default, other == printer name
{

	OSStatus		anErr, anErr2;
	AEDescList	docList;				// list of docs passed in
	long		index, itemsInList;
	Boolean		wasAlreadyOpen;
/*	
	anErr = AEGetParamDesc( inputEvent, keyDirectObject, typeAEList, &docList);
	nrequire(anErr, GetFileList);

	anErr = AECountItems( &docList, &itemsInList);			// how many files passed in
	nrequire(anErr, CountDocs);
	for (index = 1; index <= itemsInList; index++)			// handle each file passed in
	{	
		AEKeyword	keywd;
		DescType	returnedType;
		Size		actualSize;
		FSRef 		fileRef;
		FSCatalogInfo	theCatInfo;
		
		anErr = AEGetNthPtr( &docList, index, typeFSRef, &keywd, &returnedType,
						(Ptr)(&fileRef), sizeof( fileRef ), &actualSize );
		nrequire(anErr, AEGetNthPtr);

		anErr = FSGetCatalogInfo( &fileRef, kFSCatInfoFinderInfo, &theCatInfo, NULL, NULL, NULL );
/*		
		nrequire(anErr, FSGetCatalogInfo);

		if (anErr == noErr)
			anErr = DetermineWindowTypeOrOpen(&fileRef, ((FInfo*)&theCatInfo.finderInfo)->fdType, nil, nil, &wasAlreadyOpen);
			
			
*/		

//	}

	return anErr;
	
} // DoOpenOrPrint

// --------------------------------------------------------------------------------------------------------------
static pascal OSErr	DoOpenDocument(
	const AppleEvent 	*inputEvent,
	AppleEvent 	*outputEvent,
	SInt32		handlerRefCon)
{

#pragma unused (outputEvent, handlerRefCon)

	OSStatus		anErr=0;
	
    anErr = DoOpenFile(inputEvent);

	return anErr;
	
} // DoOpenDocument
void InitMac()
{
#if TARGET_API_MAC_OS8
	InitGraf(&qd.thePort);
	InitFonts();
	InitWindows();
	InitMenus();
	TEInit();
	InitDialogs(0L);
	FlushEvents(everyEvent, 0L);	
	MaxApplZone();
#endif /* TARGET_API_MAC_OS8 */
	MoreMasters();
	
	InitCursor();
#if TARGET_API_MAC_CARBON
	BitMap	screenBitMap;
	Rect	screenBits;
	Cursor theArrow;
	GetQDGlobalsScreenBits(&screenBitMap);
	screenBits = screenBitMap.bounds;
	SetCursor(GetQDGlobalsArrow(&theArrow));
#else
	Rect screenBits = qd.screenBits.bounds;
	SetCursor(&qd.arrow);
#endif /* TARGET_API_MAC_CARBON */
	
  SIOUXSettings.initializeTB = 0; // else SIOUX initialize the toolbox for us
  SIOUXSettings.toppixel = 45;
  SIOUXSettings.leftpixel = 15; 
//  SIOUXSettings.fontface = bold + italic;// or normal
  SIOUXSettings.asktosaveonclose = 0;
  
  short bas = screenBits.bottom;
  short droit = screenBits.right;  
  SIOUXSettings.columns = (short)(2.+(float)droit/8.);
  SIOUXSettings.rows = (short)(10. + (float)bas/18.);
  SIOUXSetTitle("\pfreefem+ line output"); //marche pas!!
  SIOUXSettings.fontface = normal;
  SIOUXSettings.fontid = 22;// courier;
  SIOUXSettings.fontsize = 10;
  cout << "Initmac" << endl;
  
	#define INSTALL(event, handler) \
			AEInstallEventHandler(kCoreEventClass, event, handler, 0, false)
	// AEC, changed to use the correct handler procs
//	INSTALL (kAEOpenApplication, NewAEEventHandlerUPP(DoOpenApp));
//	INSTALL (kAEReopenApplication, NewAEEventHandlerUPP(DoReopenApp));
	INSTALL (kAEOpenDocuments,   NewAEEventHandlerUPP(DoOpenDocument));
	#undef INSTALL

}
class InitilisationMac {
  static int init;
  public:
   InitilisationMac(){ InitMac();}
};

static InitilisationMac Initmac; // to call InitMac

int getprog(char* fn,int  argc, char** argvptr)
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
   else return 0; // erreur cancel
	return (2);
}
#ifdef FREEFEM

void coutmode(short i) 
{ 
   cout <<  flush;
   cerr <<  flush;
 //  if(i)(**(SIOUXTextWindow->edit)).txFace = 0;
 //  else (**(SIOUXTextWindow->edit)).txFace = 1;
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
   longjmp(environ,1);
}
void thisexit(){ myexit();}

int main (int argc, char **argv)
{
  char       *prog;
  char       fname[256];
  SIOUXSettings.sleep=1;
  argc = getprog (fname, argc, argv);
  atexit(thisexit);
  NEW_HANDLER (); // see dependent system files ({pc,x,mac}rgraph.{h,cpp})
  

  int OPTION = 0;
  if (argc == 2)
    {
        initgraphique();
        if(0==setjmp(environ))
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
   SIOUXUseWaitNextEvent = true;
   SIOUXSettings.sleep=0;
    char * wn = new  char [256];
   for (int i=0;i<256;i++)
     wn[i] = 0;
   strcpy(wn,"   -- FreeFem++ ");
   strcat(wn,StrVersionNumber());
   
   SIOUXSetTitle((unsigned char *) wn);
   int ret=15;  
                        try {                  
                          ret=mymain(argc,argv);}
                        catch( ...) { cerr << "catch exception ???";}
                        catch( Error & err) {
                          cerr  << err.what() << endl;                        
                        }

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
#if TARGET_API_MAC_CARBON
	BitMap	screenBitMap;
	Rect	screenBits;
	Cursor theArrow;
	GetQDGlobalsScreenBits(&screenBitMap);
	screenBits = screenBitMap.bounds;
	SetCursor(GetQDGlobalsArrow(&theArrow));
#else
	Rect screenBits = qd.screenBits.bounds;
	SetCursor(&qd.arrow);
#endif /* TARGET_API_MAC_CARBON */

	boundsRect.top = 45;
	boundsRect.left = (short) (15 +  (0.35 * screenBits.right));
	boundsRect.bottom = screenBits.bottom -  25;
	boundsRect.right =  screenBits.right-  25;
	if((boundsRect.bottom - boundsRect.top) < (boundsRect.right - boundsRect.left))
		boundsRect.right = boundsRect.left + boundsRect.bottom - boundsRect.top;
	else
		boundsRect.bottom = boundsRect.top + boundsRect.right - boundsRect.left;
	  grafWindow0=NewCWindow(0, &boundsRect, "\pFreeFem Graphics",true, 8, NULL, true, 0);
	ShowWindow(grafWindow0);
	BringToFront(grafWindow0);
	SelectWindow(grafWindow0);
	SetPortWindowPort(grafWindow0);
	GetPort(&grafPort0);
	height = boundsRect.bottom - boundsRect.top - 10;
	width = boundsRect.right - boundsRect.left -10;
	aspx = boundsRect.right - boundsRect.left -10;
	aspy = boundsRect.bottom - boundsRect.top - 10;
	carre = aspx == aspy;
	lacouleur = getcolor();
	CrossCurseur = GetCursor(crossCursor);
	WatchCurseur = GetCursor(watchCursor);
        GetQDGlobalsWhite(&white);
        GetQDGlobalsBlack(&black);
	
	//if( (**(wgRecord0.port.portPixMap)).pixelSize>7)
		nbcolor= 256; 
	//else 
	//	nbcolor= 2;


	ncolortable =0;
	LastColor=2;// En couleur pas defaul
	colortable=0;
	SetColorTable(2+6);
  //    TextFont(fontNum);
	TextSize(9); // small size 
	INITGRAPH = 1;
	

}


void SetColorTable(int nb)
{
  static bool greyo = !grey;
  if(!INITGRAPH) return;
   if (ncolortable == nb && greyo == grey) return;// optim
   greyo = grey;
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
     DisposeWindow(grafWindow0);
     delete [] colortable;colortable=0;
    }
  INITGRAPH=0;
}
void showgraphic()
{
  
  if (grafWindow0 != FrontWindow())
   { 
	ShowWindow(grafWindow0);
	BringToFront(grafWindow0);
	SelectWindow(grafWindow0);
	SetPortWindowPort(grafWindow0); }
	GetPort(&grafPort0);
	
}

void reffecran(void)
{
  if(!INITGRAPH) return;
    Rect rect;
    GetPortBounds(grafPort0,&rect);
    EraseRect(&rect);
 
}

int getcolor(void)
{ return 0;
}

void putpixel(int ix,int iy, int couleur)
{
  if (ncolortable>3 && couleur < ncolortable && couleur >=0 ) 
    SetCPixel(ix,iy,colortable+couleur);
}

 void plotstring(const char *s)
{  DrawText(s,0,strlen(s));
 if(psfile) fprintf(psfile,"(%s) S\n",s);
} 

int LaCouleur(){return lacouleur;}

void couleur(int c)
{ 
  if ( lacouleur == c) // small optim
    return;
  c= c > LastColor ? 1 : c; // c=Min(c,LastColor); pour noir et blanc
  lacouleur =c;
  if ( c == 0 )
    ForeColor(30);
  else if (ncolortable>3 && c < ncolortable && c >=0 ) 
    RGBForeColor(colortable+c);
  else 
   ForeColor(33);
 if (psfile)
  {
    float r=1,g=1,b=1;
    if (colortable) {
      if (c>0 && c < ncolortable)
	{
	  r =  (float) colortable[c].red /65535.F;
	  g =  (float) colortable[c].green /65535.F;
	  b =  (float) colortable[c].blue /65535.F;
	}
    }
    else if (c!=0)
      r=g=b=0;
    
    fprintf(psfile,"%.3f %.3f %.3f C\n",r,g,b);
  }
   
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
  PenSize(pepais,pepais);
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
  int xasp,yasp, getmaxx, getmaxy;
  
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
  MoveTo(newx,newy);
  if (psfile) 
   fprintf(psfile,"%d %d M\n", newx, height-newy);
  
}

void rlineto(float x, float y)
{
  int newx = scalx(x), newy = scaly(y);
  LineTo(newx,newy);
   if (psfile) 
    fprintf(psfile,"%d %d L\n", newx,height-newy);
  
}

void raffpoly(int n, float *poly)
{
  PolyHandle thePoly;
  int i;
  thePoly =OpenPoly();
   MoveTo(scalx(poly[0]),scaly( poly[1]));
    for(i=1; i<n; i++)
    LineTo(scalx(poly[2*i]),scaly( poly[2*i+1]));
    ClosePoly();
    FillPoly(thePoly,&white);
    FramePoly(thePoly);
    KillPoly(thePoly);
}
void fillpoly(int n, float *poly)
{
  PolyHandle thePoly;
  int i;
  thePoly =OpenPoly();
   MoveTo(scalx(poly[0]),scaly( poly[1]));
    for(i=1; i<n; i++)
    LineTo(scalx(poly[2*i]),scaly( poly[2*i+1]));
    ClosePoly();
    FillPoly(thePoly,&black);
    FramePoly(thePoly);
    KillPoly(thePoly);
    if (psfile) 
    {
     fprintf(psfile,"bF ");
     for (i=0;i<n;i++)
      fprintf(psfile,"%d %d ", scalx(poly[2*i]),height-scaly( poly[2*i+1]));
     fprintf(psfile,"eF\n");
    }
}

int pStrCopy (StringPtr p1, char * p2)
/* copies a pascal string `p1 into a C string */
{
	int len,i;
	
	len = (*p1++) %256;
	for(i=1;i<=len;i++) *p2++=*p1++;
	*p2 = 0;
	return 0;
}





int execute(const char* what)
{
  cout << " sorry no execute on MacOs we skip "<< what<<endl;
  system(what);
  return 1; // error
}

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
	     if (ours((WindowPtr) myEvent.message)) {
	       	BeginUpdate((WindowPtr) myEvent.message);
		    EndUpdate((WindowPtr) myEvent.message);
	     }
	    break;

}
  return char1;
}
void  viderbuff(){
   QDFlushPortBuffer(grafPort0,0);
}

char Getijc(int & x,int & y)
{   
   char char1=0;
   showgraphic();
    EventRecord		myEvent;
	int flag=1;
	HLock( (Handle) WatchCurseur);
	SetCursor(*CrossCurseur);
	HUnlock( (Handle) WatchCurseur);
	SelectWindow(grafWindow0);
   while (char1==0) {
	if (GetNextEvent(everyEvent, &myEvent) /* ,OxFFFFFFFF,h)*/) 
	  char1=HandleEvent(myEvent);
	 }
   GlobalToLocal( & myEvent.where);
   x = myEvent.where.h;
   y = myEvent.where.v;
	HLock( (Handle) WatchCurseur);
	SetCursor(*WatchCurseur);
	HUnlock( (Handle) WatchCurseur);

 //   printf("\t\t x = %d y = %d  c=%d\n", x,y,char1);
  return char1;
    
    

}

char Getxyc(float &x,float &y)
{ 
  char c;
  int i,j;
  c = Getijc( i,j);
  x = scali(i);
  y = scalj(j);
  return c;
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
   
  if(psfile==0) {printf("Erreur %s errno %d\d",filename?filename:ffff,errno);exit(1);}
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
  FontInfo 	MyFontInfo;
  GetFontInfo(&MyFontInfo);   				
 int interligne = MyFontInfo.ascent + MyFontInfo.descent + MyFontInfo.leading;
 return interligne*0.7/echy;
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

