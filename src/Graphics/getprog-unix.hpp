#include "mode_open.hpp"
#if WIN32
#include  "ff-win32.cpp"
#else
#include <unistd.h>
#endif
extern long mpirank;
extern long verbosity;
extern FILE *ThePlotStream; //  Add for new plot. FH oct 2008
// for the environ variables ...
extern const char *  prognamearg;
extern const char *  edpfilenamearg;
extern bool  waitatend;
extern bool echo_edp;
extern bool 	  NoGraphicWindow;

char * Shell_Space(const char * s);

char * Shell_Space(const char * s)
{
    const char *c=s;
    int nbspace;
    int i;
    for( i=0;i<100000;++i) 
	if (! s[i]) break;
 	else if(isspace(s[i])) ++nbspace;
    if (!(i<100000))
      {
	cerr << " Bug Shell_Space routine " <<endl;
	exit(1);
      }

#ifdef WIN32
    char * p= new char[i+1+nbspace];
    char * q=p;
    for( i=0;i<100000;++i) 
	if (! s[i]) break;
        else if(isspace(s[i])) *q++='^', *q++=s[i];
	    else *q++=s[i];
#else
    char * p= new char[i+nbspace];
    char * q=p;
    for( i=0;i<100000;++i) 
	if (! s[i]) break;
        else if(isspace(s[i])) *q++='\\', *q++=s[i];
	else *q++=s[i];
#endif
    *q++='\0';
    assert( q-p <= i+nbspace);
    return p;
}
int getprog(char* fn,int argc, char **argv)
{
  waitatend=true;  // attent 
  int ret=0;
  *fn='\0';
#ifdef WIN32
 const  int lsuffix= 4;
#else 
 const  int lsuffix= 0;
#endif
  
#ifdef PROG_FFGLUT
  const char * ffglut=PROG_FFGLUT;
#else
  const char *ffglut= "ffglut";
#endif
  const char *progffglut=0;
  const char *fileglut=0;
  bool noffglut=false;
#ifndef NODEFFFGLUT
  if(argc)
    {
      const char *prog =argv[0];
      const char *pm= strrchr(argv[0],'-');      
      if( pm )
        noffglut = ((strlen(prog)- (pm-prog)) < lsuffix+5);
      else   noffglut==  false;
      //      cout << " noffglut= " << noffglut << endl;
      //  suffix ++-glx.exe -> no ffglut
      // pm = 0= > pas de moins -> freefem++ -> ffglut
    }
#endif
  bool ch2edpdir = false;
  if(argc)
    prognamearg=argv[0];
   echo_edp=true;
  if(argc)
    for (int i=1; i<argc;i++)
      if  (ret ==0 && strcmp(argv[i],"-f")==0 && i+1 < argc  ) 
	{
	  strcpy(fn,argv[i+1]);
	  i++;	
	  edpfilenamearg=argv[i];
	  ret=1;
	}
      else if  (strcmp(argv[i],"-v")==0 && i+1 < argc) 
	{
	  verbosity = atoi(argv[i+1]);
	  i++;	
	  if(verbosity>10) printf(" verbosity : %ld\n",verbosity);
	}
      else if  (strcmp(argv[i],"-nw")==0 ) 
	{
	  noffglut=true;
	  NoGraphicWindow=true; 
	}
      else if  (strcmp(argv[i],"-ne")==0 ) // no edp 
	  echo_edp=false;
      else if  (strcmp(argv[i],"-cd")==0 ) // 
	  ch2edpdir=true;
   
      else if  (strcmp(argv[i],"-ns")==0 ) // no script  
	  echo_edp=false;
      else if  (strcmp(argv[i],"-nowait")==0 ) 
	waitatend=false;
      else if  (strcmp(argv[i],"-wait")==0 ) 
	  waitatend=true;
      else if(strcmp(argv[i],"-fglut")==0 && i+1 < argc)
	{ 
	  fileglut=argv[++i];
	  noffglut=true;
	}
      else if(strcmp(argv[i],"-glut")==0 && i+1 < argc)
	{
	  progffglut=argv[++i];
	  noffglut=true;
	  NoGraphicWindow=false;
	}
      else if(strcmp(argv[i],"-gff")==0 && i+1 < argc)
	{
	  progffglut=Shell_Space(argv[++i]);
	  noffglut=true;
	  NoGraphicWindow=false;
	}    
      else if(strcmp(argv[i],"-?")==0 )
	ret=2;
      else if( strcmp(argv[i],"-f")==0 && i+1 < argc) 
	{
	  strcpy(fn,argv[++i]);
	  ret=1;
	  edpfilenamearg=argv[i];	 
	}
      else if(ret==0)
	{
	  strcpy(fn,argv[i]);
	  edpfilenamearg=argv[i];	 
	  ret=1;
	}
if( ch2edpdir && edpfilenamearg)
  {
    int i=0;
    int l= strlen(edpfilenamearg);
#ifdef WIN32	
    const char sepdir='\\';
#else
   const char sepdir='/';
#endif    
    
    for(i=l-1;i>=0;i--)
	if(edpfilenamearg[i]==sepdir) break;
 	
    if(i>0) {
	char *dir= new char [l+1];
	strcpy(dir,edpfilenamearg);
	dir[i]=0;
	int err=0;
	if(verbosity>1) 
	    cout << " chdir '" << dir <<"'"<< endl;
#if WIN32	
	err=_chdir(dir);
#else
	err=chdir(dir);
#endif
	//cout << err << endl;
         if(err) {
	     cerr << " error : chdir  " << dir << endl;
	     exit(1);
	 }
	delete [] dir;
	
    }
      
    
  }
  if( ! progffglut && !noffglut)
    progffglut=ffglut;
  
  if(progffglut && mpirank==0)
    {
      ThePlotStream = popen(progffglut,"w");		   
      if(verbosity)
	printf(" EXEC of the plot  : %s\n",progffglut);
      if(!ThePlotStream) { cerr << "  Error popen  "<< progffglut << endl;exit(1);}
      
    }
  else if (fileglut)
    {// correction progffglut -> fileglut v3.0-2 FH.
      ThePlotStream = fopen(fileglut, MODE_WRITE_BINARY );
      if(verbosity)
	printf(" save of the plot in file : %s\n",fileglut);
      if(!ThePlotStream) 
	{
	  cerr << "  Error save file glut " << fileglut 
	       << " mode " << MODE_WRITE_BINARY<< endl;
	  exit(1);
	}
    }

#ifdef WIN32
  if(ret==0)
    {
      if ( ShowOpenDialogBox1(fn) )
	ret=1;			    
    }
#endif

  if(ret !=1) 
    {
      const char * ff = argc ? argv[0] : "FreeFem++" ;
      cout << " Syntaxe = " << ff  << " [ -v verbosity ] [ -fglut filepath ] [ -glut command ] [ -nw] [ -f] filename  [SCRIPT-arguments]\n"
	   << "        -v      verbosity : 0 -- 1000000 level of freefem output \n"
	   << "        -fglut  filepath  : the file name of save all plots (replot with ffglut command ) \n"
	   << "        -glut    command  : change  command  compatible with ffglut  \n"
           << "        -gff     command  : change  command  compatible with ffglut (with space quoting)\n"
	   << "        -nowait           : nowait at the end on window   \n"
	   << "        -wait             : wait at the end on window   \n"
	   << "        -nw               : no ffglut, ffmedit  (=> no graphics windows) \n"
	   << "        -ne               : no edp script output\n"
           << "        -cd               : Change dir to script dir\n"

	;

      if(noffglut)  cout << " without     default ffglut : " << ffglut << endl;
      else          cout << " with        default ffglut : " << ffglut << endl;
      cout   << endl;
      exit(1);
      return ret; 
    }
  if(verbosity>10) 
    cout << " file : " << fn << endl ; 

    return 1;
}
