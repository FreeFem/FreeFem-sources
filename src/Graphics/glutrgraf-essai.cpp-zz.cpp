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

#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <cassert>
#include "rgraph.hpp"

#include <pthread.h>
#include <semaphore.h>
pthread_mutex_t mutex, mutexclose;

static sem_t  k_sem_wait,k_sem_getxyc;
int k_wait_sem_wait=0;

#include "error.hpp"
using namespace std;
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


template<class T> inline T Min (const T &a,const T &b){return a < b ? a : b;}
template<class T> inline T Max (const T &a,const T & b){return a > b ? a : b;}


static  long  cube6[7][3] ={ { 65535,32000,32000},{ 65535, 65535,0},{0, 65535,0},{0, 65535, 65535},{0,0, 65535}
     , { 65535,0, 65535},{ 32000,0,0} }; 
static  long grey6[2][3] ={ {65534,65534,65534},{0,0,0} }; 

static bool grey=false;
static FILE *psfile = 0;
static FILE *psfile_save = 0;
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
static int  lacouleur,screen, width, height, currx, curry;
#define call(i) i
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

void * compile1(void * f)
{
  compile((char*) f);
  closegraphique();
  return 0;
}
/*
void * glutthread(void *argu)
{
    char ** argv =   (char **) ((void**) argu)[1];
    int &  argc =  * (int *) ((void**) argu )[0]  ;


     return return_value;     
}

void * main1(void *argu)
{
    pthread_mutex_lock(&mutex);

    char ** argv =   (char **) ((void**) argu)[1];
    int &  argc =  * (int *) ((void**) argu )[0]  ;

    if (argc <3)
     {
	cerr << " utilisation : " << argv[0] << " meshfile  solfile [rap in z ] [nbisovalue] " << endl;
	return return_value;
     }
    global=0;

    assert(argc>2);
  R rapz=1; 
  int nbiso=20;
  if (argc>3) rapz=atof(argv[3]);
  if (argc>4) nbiso=atoi(argv[4]);
  cout << " Rap z " << rapz << endl;
  Mesh Th(argv[1]);
  Rn f(Th.nv); 
  {
    ifstream fdat(argv[2]);
    assert(fdat.good());
    fdat >> f;
  } // pour ferme le fichier (la variable fdat est detruite)
    // 
  global=new Global(Th,f,100,100,rapz,nbiso);  
  pthread_mutex_unlock(&mutex);
  cout << " un lock in main " << endl;
  

   
   return return_value;
}

int main(int argc, char** argv)
{
    global=0;
    pthread_mutex_init(&mutex,NULL);
    pthread_mutex_init(&mutexclose,NULL);
    pthread_mutex_lock(&mutex);
    pthread_mutex_lock(&mutexclose);
     pthread_t tid;
    void * argu[2]={ (void *) & argc, (void*) argv};
    
    
    pthread_create(&tid,NULL,main1,(void *) argu);
    pthread_mutex_unlock(&mutex);

    glutthread(argu);
    

   
   
   void **value_ptr;
  
   pthread_join(tid,value_ptr );
  
   return 0;
}


*/

#endif
/*  
les actions :
0  init
1  fin
2  clean
3  set couleur
4  set epais
6  move
7  line
8  draw string
9  fill poly
10  draw pixel
11  videbuff
12  getxyc
*/
enum k_action { k_init,k_fin,k_clean,k_color,k_epais,k_move,k_line,k_drawstring,k_fillpoly,k_pixel,k_vide_buff,k_getxyc};
volatile int lgbuf=1024;
volatile int kbuf=0;
volatile int *buf;
int kaction[] = { 1,1,1,2,2,3,3,2,2,3,1,1 };
void k_wait()
{
    k_wait_sem_wait=1;
    sem_wait&(k_sem_wait);
}
void k_send()
{
  if(k_wait_sem_wait) sem_post(&k_sem_wait);
  k_wait_sem_wait=0;
}
void k_lock(int n)
{
  if(kbuf+n>lgbuf) {k_wait();}
  pthread_mutex_lock(&mutex);
}
void k_unlock()
{
  pthread_mutex_unlock(&mutex);
}

void addaction(int ka)
{
 assert(kaction[ka]==1 && kbuf )
  k_lock(1); 
 lgbuf[kbuf++] = ka;
  un_lock(); 
}

void addaction(int ka,int i0)
{
 assert(kaction[ka]==2 && kbuf )
  k_lock(2); 
 lgbuf[kbuf++] = ka;
 lgbuf[kbuf++] = i0;
 k_unlock(); 
}
void addaction(int ka,void *
 i0)
{
 assert(kaction[ka]==2 && kbuf )
  k_lock(2); 
 lgbuf[kbuf++] = ka;
 lgbuf[kbuf++] = i0;
 k_unlock(); 
}
void addaction(int ka,int i0,int i1)
{
 assert(kaction[ka]==3 && kbuf )
  k_lock(3); 
 lgbuf[kbuf++] = ka;
 lgbuf[kbuf++] = i0;
 lgbuf[kbuf++] = i1;
k_unlock(); 
}

void DoActions()
{

static int xo=0,yo=0;
/*  
les actions :
0  init
1  fin
2  clean
3  set couleur
4  set epais
6  move
7  line
8  draw string
9  fill poly
10  draw pixel
11  videbuff
12  getxyc
*/  pthread_mutex_lock(&mutex);
  for(int k=0,ka;k<lgbuf;)
   {
     int ka=lgbuf[k];
     int na=kaction[k];
     int i0=kaction[k+1];
     int i1=kaction[k+2];
     switch(ka) {
     case  0: // init
      break;
     case  1: // fin
      
      break;
     case  2: // clean
      		glClearColor(1.f, 1.f, 1.f, 1.0f);
		    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT); /* Clear buffer */
      break;
     case  3: // set color
      	glColor4f ((i0%256)/255.F,((i0/256)%256)/255.F,(i0/62976)/255.F,1.);
      break;
     case  4: // set epais
       glLineWidth(i0);
      break;
     case  5: // ... 
      break;
     case  6: // move
             xo=i0;yo=i1;
      break;
     case  7: // line
        glBegin(GL_LINES);
        glVertex2i(xo,yo);
        glVertex2i(i0,i1);
        xo=i0;yo=i1;
        glEnd();

      break;
     case  8: // draw string
     //  DrawCStringGL(s,fontList);
      // delete [] s;
      break;
     case  9: // fill poly
       glBegin(GL_POLYGON);
       for (int i=0;i<n;i++)
         glVertex2i(poly[2*i],poly[2*i+1]);
       glEnd();

      break;
     case  10: // draw pixel
      break;
     case  11: // vide buff
      break;
     case  12: // getxyc
      break;
     default: assert(0); // bug;
      }
     k+=na;
     
   }
  pthread_mutex_nolock(&mutex);
  k_send(); 

}
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
  width = 500;
  height =  300;
}

void closegraphique()
{
  if (INITGRAPH)
    {
      INITGRAPH = 0;
      delete [] colortable;
      closePS();
    }
  cout << " wait close " << endl;
  pthread_mutex_lock(&mutexclose);
  pthread_mutex_unlock(&mutexclose);
  cout << "  close " << endl;

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
}

void rmoveto(reel x, reel y)
{
  currx = scalx(x);
  curry = scaly(y);
}

void rlineto(reel x, reel y)
{
  int newx = scalx(x), newy = scaly(y);
  if (psfile)
    fprintf(psfile,"%d %d %d %d L\n",currx, height-curry, newx, height-newy);
  currx = newx; curry = newy;
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
{ int l = strlen(string);
 if(psfile) fprintf(psfile,"(%s) %d %d  S\n",string,currx,height-curry);
}

void showgraphic()
{
}

void x11draw3(int * ptype)
{
  int type;

  type=  *ptype;

  if (psfile) 
    switch (type) {
    case 0  : {fprintf(psfile,"[] setdash\n");break;}
    case 1  : {fprintf(psfile,"[3]  setdash\n");break;}
    default : {fprintf(psfile,"[4 1] setdash\n");break;}
    }

}  

void penthickness(int pepais)
{
  if (psfile) fprintf(psfile,"%d setlinewidth\n",pepais);
}


void x11linsrn(int * x1,int * x2,int * y1,int * y2)
  //int *x1,*x2,*y1,*y2;
{   
}

   
void viderbuff()
{
}



void cercle(reel centrex, reel centrey, reel rayon)
{
  int r = (int) (rayon * echx);
}
void reffecran()
{
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

}


int  execute (const char * str)
{ 
 cout << "exec: " << str << endl;
 return  system(str);
}

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
void rattente(int waitm)
{
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
  float s=1;
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
    fprintf(psfile," /Helvetica findfont 9 scalefont setfont\n");
    fprintf(psfile," /S {moveto show} def\n");
    fprintf(psfile," /bF  { mark} def \n");
    fprintf(psfile," /eF {newpath moveto counttomark 2 idiv {lineto} repeat closepath fill cleartomark} def\n");
    fprintf(psfile," /P { /yy exch def /xx exch def   xx xx 1 add yy yy 1 add  rec  fill } def\n");
    fprintf(psfile," 0.5 setlinewidth\n");
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
  return 5./echy;
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


void Clean() 
{
    glClearColor(1.0, 1.0, 1.0, 0.0);
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
}

static void Reshape( int lwidth, int lheight )
{   
   width=lwidth;
   height=lheight;
   // glutPostRedisplay();
}


static void Key( unsigned char key, int x, int y )
{
/*   switch (key) {
      case 27: // esc char
         exit(0);
         break;
      case '+':  
      	global->coef_dist /= 1.2;
      	break;
      case '-':  
      	global->coef_dist *= 1.2;
      	break;
      case 'g':  
      	global->theta += pi/180.;
      	break;
      case 'd':  
      	global->theta -= pi/180.;
      	break;
      case 'h':  
      	global->phi += pi/180.;
      	break;
      case 'b':  
      	global->phi -= pi/180.;
      	break;
      case 'a':
        global->dtheta = pi/180.;
        break;	
      case 's':
        global->dtheta = 0;
        break;	
      case '=':
        global->coef_dist =1;
        global->theta = 45;
        global->phi = 45;
        break;
        
      	
   }
   glutPostRedisplay();*/
}

void Display(void)
{ 
/*    Clean();

    // ALH - on doit reconstruire la liste d'affichage Žà chaque image
    // puisque les donnŽées sont susceptibles de changer.
    global->MakeListDraw();

    glutSwapBuffers();*/
}

static void Idle( void )
{

  // ALH - recalcule l'image pŽériodiquement afin d'y voir la
  // convergence de GC (ceci pourrait Žêtre optimisŽé en n'affichant une
  // nouvelle image que lorsque x change).
//  glutPostRedisplay();
}

static void Mouse( int button,int state,int x,int y )
{
 // state up or down 
 // cout << "Mouse " << button<< " " << state << " " << x << " " << y << endl;
}
static void MotionMouse(int x,int y )
{
}
void SpecialKey(int key, int x, int y)
{
}

#ifdef FREEFEM
int main (int argc, char **argv)
{
  atexit(myend);
    global=0;
    pthread_mutex_init(&mutex,NULL);
    pthread_mutex_init(&mutexclose,NULL);
    pthread_mutex_lock(&mutex);
    pthread_mutex_lock(&mutexclose);
    sem_init(&k_sem_getxyc,0,0);
    sem_init(&k_sem_wait,0,0);
    
     pthread_t tid;
    void * argu[2]={ (void *) & argc, (void*) argv};
    
    
    pthread_mutex_unlock(&mutex);
  
  int OPTION = 0;
  if (argc == 2)
    {
     glutInit(&argc , argv);
     glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH );
     height = 512;
     width = 512; 
     glutInitWindowSize(width , height);
     glutInitWindowPosition(100, 100);
     string titre = "FreeFem++ ";
     titre += argv[1] ;
     glutCreateWindow(titre.c_str());
     glutPushWindow();
     cout << "mutex lock in  glut " << endl;
     
     pthread_mutex_lock(&mutex);
     pthread_mutex_unlock(&mutex);

     glEnable(GL_DEPTH_TEST); 
     glutReshapeFunc( Reshape ); // pour changement de fenetre 
     glutKeyboardFunc( Key );    // pour les evenements clavier
     glutSpecialFunc(SpecialKey);
     glutMouseFunc(Mouse);       // pour les evenements sourie
     glutMotionFunc(MotionMouse); // les mouvements  de la sourie 
     glutDisplayFunc( Display ); // l'affichage
     glutIdleFunc( Idle );       // l'animation automatique     
     glutMainLoop();    
       printf ("PROGRAM  FreeFem 1.0 %s \n",argv[1]);
       compile (argv[1]);
       pthread_create(&tid,NULL,compile1,(void *) argv[1]);

      
    }
  else
    printf ("To launch freefem you must type freefem  and a file name\n");
  return 0;
}
#else
extern int mymain(int argc,char **argv);
int main (int argc, char **argv)
{
  int ret=mymain(argc,argv);
  return ret;
}


#endif

