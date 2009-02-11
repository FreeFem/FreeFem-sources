#include "mode_open.hpp"
#if WIN32
#include  "ff-win32.cpp"
#endif
extern long mpirank;
extern long verbosity;
extern FILE *ThePlotStream; //  Add for new plot. FH oct 2008
// for the environ variables ...
extern const char *  prognamearg;
extern const char *  edpfilenamearg;
extern bool  waitatend;
int getprog(char* fn,int argc, char **argv)
{
  waitatend=true;  // attent 
  int ret=0;
  *fn='\0';
#if WIN32
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
      // pm = 0= > pas de moin -> freefem++ -> ffglut
    }
#endif
    
  if(argc)
    prognamearg=argv[0];

  if(argc)
    for (int i=1; i<argc;i++)
      if  (ret ==0 && strcmp(argv[i],"-f")==0 && i+1 < argc  ) 
	{
	  strcpy(fn,argv[i+1]);
	  i++;	
	  ret=1;
	}
      else if  (strcmp(argv[i],"-v")==0 && i+1 < argc) 
	{
	  verbosity = atoi(argv[i+1]);
	  i++;	
	  if(verbosity>10) printf(" verbosity : %ld\n",verbosity);
	}
      else if  (strcmp(argv[i],"-nw")==0 ) 
	noffglut=true;
      else if  (strcmp(argv[i],"-nowait")==0 ) 
	waitatend=false;
      else if(strcmp(argv[i],"-fglut")==0 && i+1 < argc)
	{ 
	  fileglut=argv[++i];
	  noffglut=true;
	}
      else if(strcmp(argv[i],"-glut")==0 && i+1 < argc)
	{
	  progffglut=argv[++i];
	  noffglut=true;
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
      cout << " Syntaxe = " << ff  << " [ -v verbosity ] [ -fglut filepath ] [ -glut command ] [ -nw] [ -f] filename  \n"
	   << "        -v      verbosity :  0 -- 1000000 level of freefem output \n"
	   << "        -fglut  filepath  :  the file name of save all plots (replot with ffglut command ) \n"
	   << "        -glut    command  :  the command name of ffglut  \n"
	   << "        -nowait           : nowait at the end on window   \n"
	   << "        -nw               :  no ffglut (=> no graphics windows) \n";
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
