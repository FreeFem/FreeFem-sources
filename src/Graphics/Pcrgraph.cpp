/**************DO NOR REMOVE THIS BANNER***************/
/*  FreeFEM+ : Language for a Finite Element Method   */ 
/*  -------    Release beta:  April 1999.             */
/*  written in 4/24/2000                              */
/*  Authors: D. Bernardi, Y. Darmaillac F. Hecht,     */
/*           O. Pironneau, K.Ohtsuka                  */
/*  You may copy freely these files and use it for    */
/* teaching or research. These or part of these may   */
/* not be sold or used for a commercial purpose with- */
/* out our consent : fax (33)1 44 27 44 11            */
/* (fax)    Olivier.Pironneau@ann.jussieu.fr          */
/******************************************************/
#include "versionnumber.hpp"
#define TOSTRING1(i) #i
#define TOSTRING(i) TOSTRING1(i)

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <setjmp.h>

#include <new.h>
#include <iostream.h>
#include <fstream.h>
//#include "vect.h"
#include "error.hpp"
//#include <memory.h>

/*
** Windows includes
*/
#include <windows.h>
#include <commdlg.h>
#include <direct.h>
#include <time.h>

#include <stdio.h>
#include <fcntl.h>
#include <io.h>      //*OT  use for the console window
#include <stat.h>

#define fill thequikdrawfill

#include "rgraph.hpp"

//char *Version = "1.2.7";
void out_of_memory ();
void NEW_HANDLER (void);
void myexit(int);
void compile(char *fname);
float scali(int i);
float scalj(int j);
//int pStrCopy (StringPtr p1, StringPtr p2);
int execute(char* what);
//int DoMouseDown (int windowPart, WindowPtr whichWindow, EventRecord *myEvent);
char Getijc(int & x,int & y);
void postexit();

static int  cube6[7][3] ={ {255,0,0},{255,255,0},{0,255,0},
      {0,255,255},{0,0,255}, {255,0,255},{255,0,0} }; 
static  int grey6[2][3] ={ {255,255,255},{0,0,0} }; 
static bool grey=false;      
static int ncolortable=0;
static int LastColor=2; // LastColor=1 => Noir et Blanc >2 =>couleur

typedef struct rgb {
  BYTE r;    // red component of color
  BYTE g;  // green component of color
  BYTE b;    // blue component of color
}	rgb;

static rgb * colortable=0;
static HPEN*  hpen=0;
static HBRUSH*  hbr=0;
static  HFONT  hFont=0;
static int fontH = 0;	// The height of font
static int cstatic=1;

int getcolor();
void putpixel(int ix,int iy, int couleur);
int scalx(float x);
int scaly(float y);
void compile (char *);

void NEW_HANDLER (void){  set_new_handler (&out_of_memory);}

#define ours(w)   (w==grafWindow0)

template<class T> inline T Min (const T &a,const T &b){return a < b ? a : b;}
template<class T> inline T Max (const T &a,const T & b){return a > b ? a : b;}

char errbuf[255];
static int INITGRAPH=0;
float rayon;
static int width,height;

static FILE *psfile = 0;
static FILE *psfile_save = 0;
static float aspx, aspy, echx,echy,ech,rxmin,rxmax,rymin,rymax;
static int currx, curry;
static int carre;

static HWND        hWnd;
static WNDCLASS    rClass;
static HDC hdc;
static HANDLE hConOut=0;
const float fMinPixel = -32000; // to avoid int overflot 
const float fMaxPixel = 32000;

/* Function definitions */
BOOL Init(HINSTANCE, HINSTANCE, LPSTR, int);
int  DoMain(HINSTANCE hInstance);
LONG WINAPI OpenWindowProc1(HWND, UINT, WPARAM, LPARAM);

int getcolor();
void putpixel(int ix,int iy, int couleur);
int scalx(float x);
int scaly(float y);
void rattente (int);
BOOL inittext(VOID);

BOOL ShowOpenDialogBox(char *fileName);
BOOL CreateProjetFile(char *fileName);
char *ChangePdeToExt(char *fileName,char *ext);
BOOL mainFreeFEM();
int GetFileName(char *fullname, char  *shortname);
DWORD GetOption(char lpszCmdLine[]);
void SetConsole(HANDLE hConsole);
FILE *GetConsoleHandle(DWORD Dev);
BOOL GetConsoleBuff();
BOOL EditLog();
//void SaveMesh(Grid& t);
//void SavePlot(int D,Grid& t, Real *f);
BOOL FatalErr(char *s, int err);
//BOOL CheckSameTrig(Grid& t);

//*OT  flag for FreeFEM+/WinfFEM
#define winf_VFFEM 1
#define winf_NOWAIT 2
#define winf_NOCOLOR 4
#define winf_NOEDIT 8
#define  winf_Usage   1024

unsigned int winf_flg = 0;
// end

char FreeFemCache[256]="\0",
     shortName[256]="\0",
     fullName[256]="\0";          

void fillpoly(int n, float *poly){
   POINT *pt;
   
   pt = new POINT[n];
   for (int i=0; i < n; i++) {
      pt[i].x = scalx(poly[2*i]); pt[i].y = scaly(poly[2*i+1]);
     }
   
   if (cstatic <0  || cstatic > ncolortable)  cstatic =1;
   SelectObject(hdc,hbr[cstatic]);
   int ret = Polygon(hdc,pt,n);
   delete []  pt;
   SelectObject(hdc,hpen[n]);
    if (psfile) 
    {
     fprintf(psfile,"bF ");
     for (int i=0;i<n;i++)
      fprintf(psfile,"%d %d ", scalx(poly[2*i]),height-scaly( poly[2*i+1]));
     fprintf(psfile,"eF\n");
    }
} 

void out_of_memory ()
{
  cout << " error: operator new failed; not enough memory" << endl;
  myexit (1);
}

void erreur(char *s)
{
   ErrorExec(s,0);
 
}
void HandleWindowEvent()
{  MSG msg;
  
//if ( PeekMessage(&msg,NULL,WM_PAINT,WM_PAINT,PM_NOREMOVE)) {
if ( PeekMessage(&msg,NULL,0,0,PM_NOREMOVE)) {
    	TranslateMessage(&msg);
   		DispatchMessage(&msg);
  	}
  	
}
void SetColorTable(int n);
void raffpoly(int n, float *poly){}
int  getprog(char* fn,int argc, char** argvptr){ 
  *fn='\0';
  if ( argc>=2) strcpy(fn,argvptr[1]);
  return argc-1;}

void penthickness(int pepais){
  if (psfile) fprintf(psfile,"%d setlinewidth\n",pepais);
}
void showgraphic(){ShowWindow(hWnd, SW_SHOW ); } // UpdateWindow(hWnd);}

void thisexit(){ myexit(0);}

    
void plotstring(const char *s)
{
    static  HFONT  hOldFont;

    if (s == 0) return;
    hOldFont = (HFONT)::SelectObject(hdc,hFont);
    //::GetTextMetrics(hdc, &tm);
    ::TextOut(hdc, currx, curry-fontH, s, strlen(s));
    ::DeleteObject(hOldFont);

	  if(psfile) {
	     fprintf(psfile,"%d %d M\n", currx,height-curry);
	     fprintf(psfile,"(%s) S\n",s);
	  }
}

void rmoveto(float x, float y)
{
	currx = scalx(x);
	curry = scaly(y);
}         

void rlineto(float x, float y)
{
  HandleWindowEvent();
  int newx = scalx(x), newy = scaly(y);
  MoveToEx(hdc,currx,curry,NULL);
  LineTo(hdc,newx,newy);
  if (psfile) {
    fprintf(psfile,"%d %d M\n", currx,height-curry);  
    fprintf(psfile,"%d %d L\n", newx,height-newy);
  }
  currx = newx; curry = newy;
}

void cadre(float xmin,float xmax,float ymin,float ymax)
{
	rxmin = xmin;
	rxmax = xmax;
	rymin = ymin;
	rymax = ymax;
	echx=aspx/(xmax-xmin);
	echy=aspy/(ymax-ymin); 
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
	RECT rc;
	GetClientRect(hWnd, &rc);

	width = rc.right - rc.left;
	height = rc.bottom - rc.top;

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

int scalx(float x)
{
  return (x - rxmin) * echx;
}

int scaly(float y)
{
  return (rymax - y) * echy;
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
//  putpixel(scalx(x), scaly(y), LastColor);
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

void SetRGBpen(int n)
{
    SelectObject(hdc, hpen[n]);
}

void SetColorTable(int nb)
{
 		RECT rc;
 		GetClientRect(hWnd, &rc);
 		aspx = (float)(rc.right - rc.left);
 		aspy = (float)(rc.bottom - rc.top);
		
	if (winf_flg & winf_NOCOLOR) return;
	nb=Max(nb,23);
  	if(nb<2) nb = 2;
  	if (ncolortable == nb) return;
        if(LastColor>1) LastColor=nb-1;
  	 
 	if (hpen) for(int i=0; i<ncolortable;i++)
 					DeleteObject(hpen[i]);
 	delete [] hpen;
  	if (hbr) for(int i=0; i<ncolortable;i++)
  					DeleteObject(hbr[i]);
  	delete [] hbr;
	if (hFont) DeleteObject(hFont);
		    if(colortable) 
		       delete [] colortable;
        	colortable = new rgb[nb+1];
 		    ncolortable = nb;
        int k=0;
		    colortable[k].r = 255;
		    colortable[k].g = 255;
		    colortable[k++].b = 255;
		    colortable[k].r = 0;
		    colortable[k].g = 0;
		    colortable[k++].b = 0;
    if (nb>2) 
    { 
       	nb -= 2;
       	for (long i0=0;i0<nb;i0++,k++)
        {  
           long  i6 = i0*6;
           long  j0 = i6/nb;// in 0..6
           long  j1 = j0+1;// in 1..6
           long  k0 = i0 - (nb*j0)/6L;
           long  k1 = (nb*j1)/6L-i0;
           long  kk = (k0+k1);
           
         if (! grey)
           {
           colortable[k].r  = ((cube6[j1][0]*k0+cube6[j0][0]*k1)/kk);
           colortable[k].g  = ((cube6[j1][1]*k0+cube6[j0][1]*k1)/kk);
           colortable[k].b  = ((cube6[j1][2]*k0+cube6[j0][2]*k1)/kk);
           }
          else 
           {
           kk=nb-1;
           k1 =  i0;
           k0 = nb - i0 -1;
           colortable[k].r   = ((grey6[0][0]*k0+grey6[1][0]*k1)/kk);
           colortable[k].g   = ((grey6[0][1]*k0+grey6[1][1]*k1)/kk);
           colortable[k].b   = ((grey6[0][2]*k0+grey6[1][2]*k1)/kk);
           }
           
           /*
 		   colortable[k].r = ((cube6[j1][0]*k0+cube6[j0][0]*k1)/kk)%256;
 		   colortable[k].g = ((cube6[j1][1]*k0+cube6[j0][1]*k1)/kk)%256;
 		   colortable[k].b = ((cube6[j1][2]*k0+cube6[j0][2]*k1)/kk)%256;
 		   */
        }         
    }
     else 
      {	ncolortable  =2;}

   	hpen = new HPEN[ncolortable];
   	hbr = new HBRUSH[ncolortable];
   	
   	for(int i=0; i<ncolortable;i++)
   	{
   		hpen[i] = CreatePen(PS_INSIDEFRAME, 1,RGB(colortable[i].r, 
    					colortable[i].g, colortable[i].b));
   		hbr[i] = CreateSolidBrush(RGB(colortable[i].r, 
  						colortable[i].g, colortable[i].b));
    }
}

void couleur(int c)
{ 
  if(c!=cstatic) 
  { 
    c= c > LastColor ? 1 : c; // c=Min(c,LastColor); pour noir et blanc

    if (!(winf_flg&winf_NOCOLOR)) {
    if (c>=0 && c < ncolortable)
    	cstatic = c;
    else cstatic = 1;  
	SetRGBpen(cstatic);
  	}
  // else  SetRGBpen(1);
  }
  if (psfile)
  {
    float r=1,g=1,b=1;
    if (colortable) {
      if (c>0 && c < ncolortable)
  {
    r =  (float) colortable[c].r /255.;
    g =  (float) colortable[c].g /255.;
    b =  (float) colortable[c].b /255.;
  }
    }
    else if (c!=0)
      r=g=b=0;
    
    fprintf(psfile,"%.3f %.3f %.3f C\n",r,g,b);
  }
}
int LaCouleur(){return cstatic;}

//* Control on the graphic window
void rattente(int waitm)
{
   int i=0, j=0;
   char c;
   if (waitm)
     if(!(winf_flg&winf_NOWAIT)) c = Getijc(i,j);
}

char Getijc(int & x,int & y)
{
  char char1=' ';
  if(!INITGRAPH) 
    {
      x = 0;
      y = 0;   
      return char1;
    }
   int  cont=1;
   POINT xy;
    xy.x =0;
    xy.y =0;
  MSG msg;      
	SetWindowText(hWnd,"Click mouse to continue");
  do
   {   
      GetMessage(&msg,hWnd,0,0);// all message
      GetCursorPos(&xy);
      switch (msg.message)
       {
         case WM_LBUTTONDOWN:char1=char(251), cont=0;
         		break; 
         // with shift 248                 
         case WM_RBUTTONDOWN:char1=char(253), cont=0;
         break;  
         // with shit 250
         // if the 2 buttom, 252, et shith 249;
         case WM_CLOSE:	myexit(2);
         case WM_DESTROY:   myexit(3);
         case WM_CHAR: char1 = (TCHAR)msg.wParam; cont = 0; break;
               //case WM_KEYDOWN: char1 = (TCHAR)msg.wParam;  cont=0; break;
       default:
  					TranslateMessage(&msg);
     				DispatchMessage(&msg);
       	break;
       }
    }
    while (cont); 
//    ScreenToClient(hWnd,&xy);
 	ShowWindow(hWnd, SW_SHOW );
 // SetWindowPos(hWnd, HWND_TOPMOST, 0, 0, 0, 0, SWP_NOMOVE | SWP_NOSIZE);
    SetWindowText(hWnd,"FreeFem++ "  TOSTRING( VersionFreeFempp )  " works...");
 	RECT rc;
    ScreenToClient(hWnd,&xy);
 	GetClientRect(hWnd, &rc);
 
 	x = xy.x-rc.left;
 	y = xy.y-rc.top;
 //	cout << " x = " << x << " y = " << y  << " char = " << ((unsigned char)char1 > 127 ? '*': char1) << ")" << endl;
 	return char1;
}

char Getxyc(float &x,float &y)
{ 
  char c=' ';
  int i=0,j=0;
  if(!(winf_flg&winf_NOWAIT)) c = Getijc( i,j);
  x = scali(i);
  y = scalj(j);
  //rattente(1);
  return c;
}

//* clear the screen with white
void reffecran(void)
{     
 HBRUSH hbr;
 RECT rc;

 GetClientRect(hWnd, &rc);
 hbr = CreateSolidBrush(RGB(255, 255, 255));
 FillRect(hdc,&rc,hbr);
 DeleteObject(hbr);
}

BOOL ShowOpenDialogBox(char *fileName)
{
  OPENFILENAME ofn; 
  char szDirName[256];   
  char *strFilter="PCgFEM Files (*.edp)\0*.edp\0All Files (*.*)\0*.*\0\0"; 
  
  memset(&ofn, 0, sizeof(OPENFILENAME));
  getcwd(szDirName,sizeof(szDirName));
  ofn.lStructSize = sizeof(OPENFILENAME);
  ofn.hwndOwner = NULL;
  ofn.lpstrFilter = strFilter;
  ofn.lpstrFileTitle = fileName;
  ofn.nMaxFileTitle = 80;
  ofn.lpstrInitialDir=szDirName;
  ofn.lpstrTitle ="Choose you freefem '*.edp' File";
  ofn.Flags=OFN_SHOWHELP|OFN_PATHMUSTEXIST|OFN_FILEMUSTEXIST;
  
  return GetOpenFileName(&ofn);
} 

void coutmode(short r) { ;}// will be done later
   

void initgraphique(void)       
{ 
  if (INITGRAPH) return;
  hdc=GetDC(hWnd);
  hpen=0;
  SetColorTable(2+6);

  RECT rc;
  GetClientRect(hWnd, &rc);
  aspx = (float)(rc.right - rc.left);
  aspy = (float)(rc.bottom - rc.top);
  width = rc.right - rc.left;
  height = rc.bottom - rc.top;
  carre = aspx == aspy;
  // Define the font style
    LOGFONT lf;
    TEXTMETRIC tm;
		HFONT hFont, hOldFont;

    memset(&lf, 0, sizeof lf);
    lf.lfHeight = -9;
    lstrcpy(lf.lfFaceName,"Arial");
    lf.lfOutPrecision = OUT_TT_PRECIS;
    lf.lfClipPrecision = CLIP_DEFAULT_PRECIS;
    lf.lfQuality = PROOF_QUALITY;
    lf.lfPitchAndFamily = FF_SWISS | VARIABLE_PITCH;
    hFont = ::CreateFontIndirect(&lf);

    hOldFont = (HFONT)::SelectObject(hdc,hFont);
    ::GetTextMetrics(hdc, &tm);
    ::DeleteObject(hOldFont);
    fontH =  (tm.tmHeight + tm.tmExternalLeading)*0.6;
  // end of font style
  INITGRAPH = 1;
  // cout << flush << "end  inigraphique " << endl;
}

void closegraphique(void)
{ 
	if(INITGRAPH) {
		if(hpen) 
		   DeleteObject(hpen), delete [] colortable;
		if (hbr)   DeleteObject(hbr), 
    INITGRAPH =0; // before DestroyWindow to avoid loop 
    ReleaseDC(hWnd,hdc); 
//    DestroyWindow(hWnd);
	}
}

void openPS(const char *filename )
{ 
  RECT rc;
  GetClientRect(hWnd, &rc);
  width = rc.right - rc.left;
  height = rc.bottom - rc.top;

  closePS();
  time_t t_loc;
  float s=0.5;
    const int shiftx=50,shifty=50;

  time(&t_loc);
  printf(" Save Postscript in file '%s'\n",filename?filename:"freefem.ps"),
  psfile=fopen(filename?filename:"freefem.ps","w");
  if(psfile==0) {printf("Erreur %s \d",filename);exit(1);}
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
    fprintf(psfile_save,"showpage\n");//fprintf(psfile,"showpage\n");
    fclose(psfile_save);//fclose(psfile);
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
     else   psfile=0;
   }

// Various works when the program will end
void myexit(int err)
{
 	time_t ltime;         // write the time stump in console
 	struct tm *now;
 	time(&ltime);           // write the end time
 	now = localtime(&ltime);
 	cout << "\nEnd Time: " << asctime(now) << endl;

  if (err==0) {  // normal end
		cout << "end No Error " << endl << flush ;
	}
  else  
  	cout << "end by Error (no.=" << err << ')' << endl;
  rattente(1);
	
  if (GetConsoleBuff()==FALSE)
 		FatalErr("Log file creat err !",0);
 	if (!(winf_flg&winf_NOEDIT)) EditLog();

  if (INITGRAPH)
   closegraphique();
 
 	FreeConsole();
	PostQuitMessage(0);
 	exit(err);
}

// initialize the console
BOOL inittext(VOID)
{
  FILE *fp;
  OSVERSIONINFO osVer; // for GetVersionEx()

  osVer.dwOSVersionInfoSize = sizeof(osVer);
  GetVersionEx(&osVer);
  if (osVer.dwPlatformId == VER_PLATFORM_WIN32s) {
    MessageBox(NULL, 
        "This FreeFEM++ cannot run on Windows 3.1.\n"
        "This application will now terminate.",
        "Error: Windows NT or Windows 95 Required to Run",  MB_OK );
        return FALSE;       // Console API is not able in Windows 3.1 
    }

   FreeConsole();          // If the console is already used
   AllocConsole();         // Use the console API

   // get the standard output
   if((fp = GetConsoleHandle(STD_OUTPUT_HANDLE)) == NULL){
      FreeConsole();
      return FALSE;
   }
   *stdout = *fp;

   ios::sync_with_stdio();

   SetConsoleTitle("FreeFEM++ v" TOSTRING( VersionFreeFempp )  " console");
   return TRUE;
}

//*------- Modules for MS-Windows
//*OT  12/3/1999
//* Get the buffer of the console
//* The buffer is stored in the filename.log 
BOOL GetConsoleBuff()
{
  CONSOLE_SCREEN_BUFFER_INFO csbi; //* to get buffer info
  GetConsoleScreenBufferInfo(hConOut, &csbi);

  COORD coordLine = {0,0};
  CHAR *szLine;  //* buffer to read from the console (a line)
  DWORD dwCharsRead;
  char fname[255];
  FILE *fp;
  strcpy(fname,ChangePdeToExt(shortName,"log"));
  if ((fp = fopen(fname,"w"))==NULL) {
		perror(fname);
     return FALSE;
	}
	
  szLine = (CHAR *)malloc((csbi.dwSize.X+1) * sizeof(CHAR));
  for (int i=0; i<csbi.dwCursorPosition.Y; i++) {
  	if (ReadConsoleOutputCharacter(hConOut, szLine,
              csbi.dwSize.X, coordLine, &dwCharsRead)== FALSE) {
       				perror("ReadConsoleOutputCharacter");
              return FALSE;
    }
    int j=csbi.dwSize.X-1;
    while ((szLine[j] == ' ') && (j > 0)) szLine[j--] =0;
    if (j < csbi.dwSize.X-1) szLine[j+1] = '\n'; 
    fprintf(fp,"%s",szLine);
    coordLine.Y++;
  }
  fclose(fp);
  return TRUE;
}

//*OT  12/3/1999
//* Open the filename.log by the editor
//* default editor is notepad.exe
//* Using variable "ffemEd", we can change the editor 
BOOL EditLog()
{
  char *editor, fname[256], cmdLine[255];
  
  strcpy(fname,ChangePdeToExt(shortName,"log"));
  editor = getenv("ffed");
  if (editor == 0)
  	sprintf(cmdLine,"notepad.exe %s",fname);
	else
	  sprintf(cmdLine,"%s %s",editor,fname);

  if (WinExec(cmdLine,SW_SHOWNORMAL) < 31) {
  	sprintf(errbuf,"Cannot execute [%s]",cmdLine);
  	FatalErr(errbuf,99);
  	return FALSE;
  }
  FreeConsole();
  return true;
}

void Usage()
{
   cout << "Usage: freefem++ [options]" << endl;
   cout << "Select a program file by the dialog box if option is omitted.\n[option]" << endl;
   cout << "-f filename: Run the program file \"filename\"." << endl;
   cout << "    In this mode, all plotted datas are stored in the \".\\cache\"." << endl;
   cout << "    The stored datas are used in \"WinfFEM\" (IDE for freefem+)." << endl;
   cout << "    You can get this from  <http://barnard.cs.hkg.ac.jp>." << endl;
   cout << "-s    : No wait at end." << endl;
   cout << "-b    : Do not use the color" << endl;
   cout << "-n    : Do not open the log file at end. The editor is the notepad if you do not" << endl;
   cout << "        set \"ffed=[name of editor]\" in environments." << endl;
   cout << "-h    : Display the usage (this)." << endl;
}

// freefem+  arg1  arg2 arg3
// Hack the args and analysis 
int StoreFname(char Line[], int len)
{
   	int i;
   	char msg[256]; char *ext;
   	for (i=0; i<len; i++)
      if (Line[i] != ' ') fullName[i] = Line[i];
      else break;
	fullName[i] = '\0';
 	
 	ofstream check(fullName,ios::in);
 	if (!check.is_open()) {
 	   sprintf(msg,"%s does not exist!",fullName);
 	   FatalErr(msg,-1);
 	}
 	else check.close();
 	   
	ext = strrchr(fullName,'.'); ext++;
	if (toupper(*ext) != 'E' || toupper(*(ext+1)) != 'D' || toupper(*(ext+2)) != 'P') {
 	   sprintf(msg,"%s is not program!",fullName);
 	   FatalErr(msg,-1);
 	}
	
   	GetFileName(fullName,shortName);
    return i;
}

// freefem+  arg1  arg2 arg3
// Hack the args and analysis 
DWORD GetOption(char lpszCmdLine[])
{
  int i = 0;
  int CmdLen = strlen(lpszCmdLine);
  DWORD dwStyle = WS_OVERLAPPEDWINDOW;

  while (i < CmdLen) { 
    while (lpszCmdLine[i] == ' ')
              i++;
  	if (lpszCmdLine[i] == '-') {
    	i++;
    	switch(lpszCmdLine[i]) {
    	case 'f':
      	i++;
      	while (lpszCmdLine[i] == ' ')
              ++i;
      	i += StoreFname(&lpszCmdLine[i],CmdLen-i);
  		winf_flg |= winf_VFFEM;
      	break;
    	case 's':  // not wait at end of execution
      		winf_flg |= winf_NOWAIT; ++i;
      		break;
    	case 'b':	// no color
      		winf_flg |= winf_NOCOLOR; ++i;
      		break;
    	case 'n':
      		winf_flg |= winf_NOEDIT; ++i;
      		break;
    	case 'h':
      		winf_flg |= winf_Usage; ++i;
      		break;
    	default:
      		while (lpszCmdLine[i]!=' ' && (i < CmdLen))
        	i++;
    	}
  	}
  	else  {
  		i += StoreFname(&lpszCmdLine[i],CmdLen-i);
  		break;
  	}
  }
  return 0;
}

/*
 * Init     
 *     Initialization for the program is done here:
 *     1)  Register the window class (if this is the first instance)
 *     2)  Create the desktop window for the app.
 *     3)  Show the desktop window in the manner requested by the User.
 *
 */
BOOL Init(HINSTANCE hInstance,   HINSTANCE hPrevInstance,
    LPSTR  lpszCmdLine, int    nCmdShow) 
{
  DWORD dwStyle = WS_OVERLAPPEDWINDOW;

  if (!hPrevInstance)
  {
    /*  Register Class for First Overlapped Window  */
    rClass.lpszClassName = "FreeFem++" ;
    rClass.hInstance     = hInstance;
    rClass.lpfnWndProc   = OpenWindowProc1;
    rClass.hCursor       = LoadCursor(NULL, IDC_ARROW);
    rClass.hIcon         = LoadIcon(NULL, IDI_APPLICATION);
    rClass.hbrBackground = (HBRUSH)GetStockObject(WHITE_BRUSH);
    rClass.style         = 0x4000;
    rClass.cbClsExtra    = 0;
    rClass.cbWndExtra    = 0;

    if (!RegisterClass( &rClass))
      return FALSE;
    }
  int dd=600;
  int ddx0=200;
  int ddy0=30;
//  long dwFlags;
/*  DEVMODE dev_mode = {0};
  if(!EnumDisplaySettings(NULL,ENUM_CURRENT_SETTINGS,&dev_mode))
	{
	   cout  << " screen size ??  " << dev_mode.dmPelsWidth << " x " << dev_mode.dmPelsHeight << endl;
	   dd = Min(dev_mode.dmPelsWidth*0.7,dev_mode.dmPelsHeight*0.9);
	   ddx0 = dev_mode.dmPelsWidth*0.28;
	   ddy0=dev_mode.dmPelsHeight*0.05;	   
	}
  else cout << " Error EnumDisplaySettings => no screen size " << endl;
  
  */
  int sx = GetSystemMetrics(SM_CXSCREEN);
  int sy = GetSystemMetrics(SM_CYSCREEN);
  dd = Min(sx*0.7,sy*0.9);
	   ddx0 = sx*0.28;
	   ddy0=sy*0.05;	   
  
  //cout << " Screen Size " << sx << " x " << sy << endl;
 // Rectangle  ss=Get_VirtualScreen();
 //  dd=(Abs(ss.get_Top-ss.get_Bottom())*90)/100;
 
  GetOption(lpszCmdLine);
  hWnd = CreateWindow("FreeFEM++",
      "FreeFEM++ v"  TOSTRING( VersionFreeFempp ) " for Windows",
      dwStyle,ddx0,ddy0,dd,dd,/*
      CW_USEDEFAULT,
      CW_USEDEFAULT,
      CW_USEDEFAULT,
      CW_USEDEFAULT,*/
      NULL,
      NULL,
      hInstance,
      NULL);

  
  if (*fullName == '\0' && (winf_flg != winf_Usage)) { // in command line, there is no filename
    if (ShowOpenDialogBox(shortName)==FALSE) {
       exit(0);
    }
    strcpy(fullName,shortName);
  }
  
  if (inittext()==FALSE) 
    myexit(1);
  else if (winf_flg & winf_VFFEM ) { // create only cache, option "-f" is given
   if (!_getcwd(FreeFemCache,MAX_PATH)) {
     FatalErr("Fail to get current path",-1);
   }
   strcat(FreeFemCache,"\\cache\\");

   if (_chdir(FreeFemCache)) {  // check the cache directory
     if (mkdir(FreeFemCache,0)) {
       sprintf(errbuf,"Fail to create the directory %s",FreeFemCache); 
       FatalErr(errbuf,-1);
     }
   }
   else (_chdir("..\\"));  // already created
   return TRUE;
  };
  
  return TRUE;
}

/* OpenWindowProc1 - Handles messages for the main window.
 *     Parameters:
 *         hWnd    - Handle to Window which message is delivered to.
 *         msgID   - ID number of message
 *         wParam  - 16-bit parameter
 *         lParam  - 32-bit parameter
 *
 */
LONG WINAPI OpenWindowProc1(
    HWND    hWnd,
    UINT    wMsgID,
    WPARAM  wParam,
    LPARAM  lParam)
{
  switch (wMsgID) {
    case WM_DESTROY:
      PostQuitMessage(0);
 			DestroyWindow(hWnd);
      break;
    default:
      return DefWindowProc(hWnd, wMsgID, wParam, lParam);
  }

  return 0;
}


//*OT  29/12/98
// Routines and functions for WinfFEM
int  chkCacheDir();
BOOL TestProjetPresence(char *shortName);
BOOL CreateProjetFile(char *shortName);
BOOL SaveLogFile(char *fileName);
void GetOption(int argc, char *argv[]);
FILE *projet=NULL;
// end 


int WINAPI WinMain(HINSTANCE  hInstance,
        HINSTANCE hPrevInstance,
        LPSTR  lpszCmdLine,  // int argc, char *argv[]
        int    nCmdShow)
{
  MSG msg;
  
  LPTSTR cmd = GetCommandLine();
  if (Init(hInstance, hPrevInstance,lpszCmdLine,nCmdShow)) {
		// main after checking options
		if (mainFreeFEM() == FALSE)
		   myexit(99);	// exit with error
		else myexit(0);

   	while (GetMessage(&msg,NULL,0,0)) {
    	TranslateMessage(&msg);
   		DispatchMessage(&msg);
  	}
  	myexit(msg.wParam);  // exit(msg.wParam);
	}
	return -1;
}
// the real main 


extern int mymain(int argc,char **argv);

// main() in FreeFEM+ for PCs
BOOL mainFreeFEM()
{
 char prjName[256]; 

 if (winf_flg & winf_VFFEM) {	// given by "-f filename"
  strcpy(prjName,ChangePdeToExt(shortName,"prj"));

  if (strcmp(FreeFemCache,"")!=0)
      if (CreateProjetFile(prjName)==FALSE)
        FatalErr(prjName,-1);
 }   
 
 cout << "Welcome to freefem++ v " << TOSTRING( VersionFreeFempp )  <<endl;
 cout << "Program file [" << fullName <<']'<< endl;
 time_t ltime;         // write the time stump in console
 struct tm *now;
 time(&ltime);
 now = localtime(&ltime);
 cout << "Start Time: " << asctime(now) << endl;
 ShowWindow(hWnd,SW_SHOW);
 initgraphique();
 int argc;
 char * argv [3];
 argc = 2;
 argv[0]= "FreeFem++";
 argv[1]= fullName;
// int main (int  argc, char **argv)
 int ret=1;
 try {
 	//compile(fullName);
 	ret=mymain(argc,argv);
 	cout << fullName << endl;
 	rattente(1);
 	
 	rattente(1);
 	} catch(Error &e){cout<<" error "<<e.what(); myexit(1); };                   
 
 if (projet!=NULL)
  fclose(projet);

 SetWindowText(hWnd,"End of FreeFEM++");
 if (winf_flg & winf_NOWAIT) {	// option "-s" 
   	myexit(0);
 };
  return true;
}

// error without console
BOOL FatalErr(char *s, int err)
{
  int ret;
  UINT Style = MB_SYSTEMMODAL|MB_ICONEXCLAMATION;

  if (err==0) Style |= MB_YESNO;
  else Style |= MB_OK;

  ret = MessageBox(NULL,s,"Information from FreeFEM++",Style);

  if (!err & (ret ==IDYES)) return TRUE;
  else if (!err & (ret == IDNO)) return FALSE;
  else if (err != 99) exit(err);
  return TRUE;
}


//*OT
// Modules for debug console and log
void SetConsole(HANDLE hConsole)
{
	SMALL_RECT srctWindowRect;
	COORD  coordScreen;

	srctWindowRect.Left = 10;
	srctWindowRect.Top = 10;   
	srctWindowRect.Bottom = 500;
	srctWindowRect.Right = 320;
	SetConsoleWindowInfo(hConsole, TRUE, &srctWindowRect);
	coordScreen.X = 80;		// buffer size 80x1000
	coordScreen.Y = 1000;
	SetConsoleScreenBufferSize(hConsole,coordScreen);
}

// the console handle to stdout 
FILE *GetConsoleHandle(DWORD Device)
{
  int Crt;
  FILE *Console = NULL;

  if ((hConOut = GetStdHandle(Device)) != INVALID_HANDLE_VALUE) {
    Crt = _open_osfhandle((long)hConOut, _O_TEXT);
    if (Device == STD_INPUT_HANDLE) 
       Console = _fdopen(Crt, "r");
    else Console = _fdopen(Crt, "w");
    setvbuf(Console, NULL, _IONBF, 0);
    SetConsole(hConOut);
  }

  return Console;
}


int GetFileName( char *fullname,
    char  *shortname)  // filename
{
    int   i, j, k;
    int   tail=0;

    ifstream test(fullname,ios::in);
    if (!test.is_open()) {
      cout << "File " << fullname << "do not exist !" << endl;
         return FALSE;
    }

    strcpy( shortname , "\0" )  ;

    tail = strlen( fullname ) ;

    if  ( tail == 0 )   return -1 ; // return by nothing

    for ( i = tail - 1 ; i >= 0 ; i-- ) { // loop 1
       if  ( fullname[i] == '\\' ) {
       for ( j = i+1, k=0 ; j < tail ; j++, k++ )
          *(shortname + k) = *(fullname + j ) ;
         *(shortname + k) = '\0';
       break;
      }
  }
  if (i == -1)
    strcpy(shortname,fullname);

   return 0  ; // OK!
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

char *ChangePdeToExt(char *fileName,char *ext)
{
 int len;
  
 len=strlen(fileName);        
 char *file = new char[len+1];
 for(int i=0; i<len; i++) *(file+i) = *(fileName+i);
 file[len-4]='.';
 file[len-3]=ext[0];
 file[len-2]=ext[1];
 file[len-1]=ext[2];
 file[len]='\0';
 return file;
} 

BOOL CreateProjetFile(char *fileName)
{
 char chemin[256];

 strcpy(chemin,FreeFemCache);
 strcat(chemin,ChangePdeToExt(shortName,"prj"));
 if (projet=fopen(chemin,"w"),!projet)
     return FALSE;
 fprintf(projet,"FFF@WinfFEM@FFF\n");                         
 return TRUE;
}

//* the module for graph.cpp
int NbMeshTotal=1, NbPlotTotal=1;

float GetHeigthFont()
{
    return (float)fontH/echy;
}

//* store all data created in graph.cpp

// Chack the mesh which will be stored in the cache
// return TRUE  (if it is same)
// else return FALSE  
/*
BOOL CheckSameTrig(Grid& t)
{
  static struct {  // store the information of the mesh
    float x;  float y;  int w;
  }  p[3] = { {0,0,0}, {0,0,0}, {0,0,0}};  // three virteces
  static int nv = 0, nt = 0;  // numbers of the virteces and the triangles

  if ((nv != t.nv) || (nt != t.nt))
    goto SET;
  if ((t.v[0].x != p[0].x) || (t.v[0].y != p[0].y) || (t.v[0].where != p[0].w))
    goto SET;
  if ((t.v[nv/2].x != p[1].x) || (t.v[nv/2].y != p[1].y) || (t.v[nv/2].where != p[1].w))
    goto SET;
  if ((t.v[nv-1].x != p[2].x) || (t.v[nv-1].y != p[2].y) || (t.v[nv-1].where != p[2].w))
    goto SET;
  else return TRUE;

SET:
 nv = t.nv;  nt = t.nt;
 p[0].x = t.v[0].x; p[0].y = t.v[0].y; p[0].w = t.v[0].where; 
 p[1].x = t.v[nv/2].x; p[1].y = t.v[nv/2].y; p[1].w = t.v[nv/2].where; 
 p[2].x = t.v[nv-1].x; p[2].y = t.v[nv-1].y; p[2].w = t.v[nv-1].where; 
 return FALSE;
}


// ohtsuka 8/23/98
 void SaveMesh(Grid& t)
{
 char chemin[256],meshName[256];
 int i,j=0;
 
 if (!(winf_flg&winf_VFFEM)) return;
 if (CheckSameTrig(t)==TRUE)
   return;

 strcpy(chemin,FreeFemCache);
 sprintf(meshName,"%d-%s",NbMeshTotal++,ChangePdeToExt(shortName,"msh"));
 strcat(chemin,meshName);
 
 ofstream mesh(chemin,ios::out);
 if (mesh.is_open())  {
   mesh << "FFF@WinfFEM_MESH@FFF" << endl;
   fprintf(projet,"%s\n",meshName);
   mesh << t.nv << "	" << t.nt << endl;
   for( i=0; i<t.nv; i++ ) {
      mesh << t.v[i].x <<"	"<< t.v[i].y <<"	"
           << t.v[i].where << endl;
   }

   for( i=0; i<t.nt; i++ ) {	
       mesh << t.no(t.t[i].v[0]) <<"	"<< t.no(t.t[i].v[1]) <<"	"<< t.no(t.t[i].v[2]) <<"	"<< j<<endl;
   }
 } else {
   cout << "Unable to SAVE MESH for WinfFEM !" << endl;
   return;
  }
  NbPlotTotal=1;
  mesh.close();
} 

void SavePlot(int D,Grid& t, Real *f)
{                                    
 char chemin[256],plotName[256];
 int i;
 
 if (!(winf_flg&winf_VFFEM)) return;
 if (f == NULL) return;
 strcpy(chemin,FreeFemCache);
 sprintf(plotName,"%d-%d%s",NbMeshTotal-1,NbPlotTotal,ChangePdeToExt(shortName,"fnc"));
 strcat(chemin,plotName);

 ofstream plot(chemin,ios::out);
 if (plot.is_open()) {
   fprintf(projet,"%s\n",plotName);
   plot << "FFF@WinfFEM_PLOT@FFF" << endl;          
   plot << NbMeshTotal-1 << "-" << ChangePdeToExt(shortName,"msh") << endl;
   plot << D << endl;
   plot << t.nv << endl;
   for (i=0; i<t.nv; i++) 
     plot << t.v[i].where << "  " << f[i] << endl;
 } else {
   cerr << "Unable to SAVE PLOT for WinfFEM !" << endl;
   return;
 }
 NbPlotTotal++;
}                         
*/
void  viderbuff(){;}

//*OT July 7 2000
//This module is used in analyse.cpp "system"
//
#include <shellapi.h>
#define MAXPATH  256
char * getOp(const char *what)
{
    int   tail=0, len=0;
    char  *p;

    p = strrchr(what, '\\');
    if (p == NULL) p = (char *)what;
    while (*p && (*p != ' '))  p++;
    if (*p) *p++ = '\0';
    else return NULL;
    while (*p && (*p == ' '))  p++;
    if (*p) return p;
    else return NULL;
}

int execute(const char* what)
{
	char szBuffer[MAXPATH + 1];
    char *option;
    int r=0; 
   	option = getOp(what);

	if (*what) {
		STARTUPINFO si;
		PROCESS_INFORMATION pi;
		ZeroMemory( &si, sizeof(STARTUPINFO) );
		ZeroMemory( &pi, sizeof(PROCESS_INFORMATION) );
		si.cb=sizeof( STARTUPINFO );
		si.dwFlags = STARTF_USESHOWWINDOW;
		si.wShowWindow = SW_SHOWNORMAL;

		CreateProcess(what,option,NULL,NULL,FALSE,CREATE_NEW_CONSOLE,NULL,NULL,&si,&pi );
		if( pi.hProcess )	{
			WaitForInputIdle( GetCurrentProcess(), INFINITE );
			DWORD dwExitCode = STILL_ACTIVE;
			while(dwExitCode == STILL_ACTIVE)	{
				WaitForSingleObject( pi.hProcess, 1000 );
				GetExitCodeProcess( pi.hProcess, &dwExitCode );
			}
			CloseHandle(pi.hProcess);
			CloseHandle(pi.hThread);
		}
		else {
	    	sprintf(szBuffer,"%s: cannot execute",what);
	    	if (option != NULL) sprintf(szBuffer,"%s  with option %s!",szBuffer, option); 
			MessageBox(NULL, szBuffer, "Error in FreeFem++", MB_OK | MB_ICONINFORMATION);
			r=1; 
		}
	}
	else r=2,MessageBox(NULL, "Error in system()", "Error in FreeFem++", MB_OK | MB_ICONINFORMATION);
  return r;
}
  void setgrey(bool gg ){grey=gg;}
  int getgrey(){ return grey;}
