#include <direct.h>
#include <string.h>
#include <cstdlib>
#include <string>
#include <cstring>
#include <iostream>
using namespace std;


const char C='"';


const char SLACH='/';
const char BACKSLACH='\\';
#ifdef PURE_WIN32
const  char dirsep=BACKSLACH, dirnsep=SLACH;
#else
const  char dirnsep=BACKSLACH, dirsep=SLACH;
#endif
string DirName(const char * f)
{
    const char *c= strrchr(f,dirsep);
    if(!c) return string("");
    else return string(f,strlen(f)-strlen(c))+dirsep;
}
int main(int argc,const char **argv)
{
  int debug=0; 
  char *dirff=0;
  const char *pp=0;
  LPWSTR buffer; //or wchar_t * buffer;
  int lbuffer = GetModuleFileName(NULL, buffer, MAX_PATH) ;
  string dirff;
  if lbuffer>0 dirff = DirName(buffer);
    
  string cmd=dirff+"freefem++.exe ";
  if(argc <=1)
  {
        cerr << " Sorry no file name "<< endl;
        cerr << " Drag and Drop the file icon on the application  icon or double clip on script file" << endl;
      cmd += " -wait -log";
      int ret= system(cmd.c_str());
      return 0; 
  }
  
  for(int i=1;i<argc;++i)
    {	
	if(strcmp("++d",argv[i])==0) 
	  debug=1;
	else {
    cmd += C;
    cmd += argv[i];
    if(!pp &&strlen(argv[i])>2) 
	if( argv[i][1]==':')
          pp= argv[i];
    cmd += C;
    cmd += " ";
	if( debug) cout << "  ffl: arg " << i << argv[i] << endl;
    }}
  if(pp)
   {
   	if( debug ) cout << "  ffl: file:" << pp << endl;  
    int i=0;
    int l= strlen(pp);
     for(i=l-1;i>=0;i--)
       if(pp[i]=='\\') break;
     dir= new char [l+1];
     strcpy(dir,pp);
     dir[i]=0;
	 if(debug) 
     cout << "  ffl:  chdir to " << dir << endl;
     _chdir(dir);
     delete [] dir;
   }
   cmd += " -wait -log";
   if(debug) 
   cout << "exec " << cmd << endl;
   int ret= system(cmd.c_str());
   return ret;
}
