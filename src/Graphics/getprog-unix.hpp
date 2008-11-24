#include "mode_open.hpp"

#ifdef WIN32

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

#endif

extern long verbosity;
extern FILE *ThePlotStream; //  Add for new plot. FH oct 2008
int getprog(char* fn,int argc, char **argv)
{
  
  int ret=0;
  *fn='\0';
  
#ifdef PROG_FFGLUT
  const char * ffglut=PROG_FFGLUT;
#elifdef WIN32
  const char * nodefffglut="nw.exe";
#else
  const char *ffglut="ffglut";
#endif
#ifdef WIN32
  const char * nodefffglut="nw.exe";
#else
  const char * nodefffglut="nw";
#endif
  const char *progffglut=0;
  const char *fileglut=0;
  bool noffglut=false;
#ifndef NODEFFFGLUT
  if(argc)
    noffglut= ( strcmp(argv[0]+strlen(argv[0])-strlen(nodefffglut),nodefffglut)==0);
#endif
    
  
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
      else if(strcmp(argv[i],"-fglut")==0 && i+1 < argc)
	fileglut=argv[++i];
      else if(strcmp(argv[i],"-glut")==0 && i+1 < argc)
	progffglut=argv[++i];
      else if(strcmp(argv[i],"-?")==0 )
	ret=2;
      else if( strcmp(argv[i],"-f")==0 && i+1 < argc) 
	{
	  strcpy(fn,argv[++i]);
	  ret=1;
	}
      else if(ret==0)
	{
	  strcpy(fn,argv[i]);
	    ret=1;
	}
#ifdef WIN32
  if(ret==0)
    {
      if ( ShowOpenDialogBox(fn) )
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
	   << "        -nw               :  no ffglut (=> no graphics windows) \n";
      if(noffglut)  cout << " without     default ffglut : " << ffglut << endl;
      else          cout << " with        default ffglut : " << ffglut << endl;
      cout   << endl;
      exit(1);
      return ret; 
    }
  if(verbosity>10) 
    cout << " file : " << fn << endl ; 

  if( ! progffglut && !noffglut)
    progffglut=ffglut;
  
  if(progffglut)
    {
      ThePlotStream = popen(progffglut,MODE_WRITE_BINARY);		   
      if(verbosity)
	printf(" EXEC of the plot  : %s\n",progffglut);
      if(!ThePlotStream) { cerr << "  Error popen  "<< progffglut << endl;exit(1);}
      
    }
  else if (fileglut)
    {
      ThePlotStream = fopen(progffglut, MODE_WRITE_BINARY );
      if(verbosity)
	printf(" save of the plot in file : %s\n",fileglut);
      if(!ThePlotStream) { cerr << "  Error save file glut " << fileglut << endl;exit(1);}
    }
    return 1;
}
