#include <direct.h>
#include <string.h>
#include <cstdlib>
#include <string>
#include <cstring>
#include <iostream>
using namespace std;


const char C='"';

int main(int argc,const char **argv)
{
  char *dir=0;
  const char *pp=0; 
  string cmd="freefem++.exe ";
  for(int i=1;i<argc;++i)
    {	
    cmd += C;
    cmd += argv[i];
    if(!pp &&strlen(argv[i])>2) 
	if( argv[i][1]==':')
          pp= argv[i];
    cmd += C;
    cmd += " ";
    }
  if(pp)
   {
    int i=0;
    int l= strlen(pp);
     for(i=l-1;i>=0;i--)
       if(pp[i]=='\\') break;
     dir= new char [l+1];
     strcpy(dir,pp);
     dir[i]=0;
     //cout << " chdir to " << dir << endl;
     _chdir(dir);
     delete [] dir;
   }
   cmd += " -wait";
   //cout << "exec " << cmd << endl;
   int ret= system(cmd.c_str());
   return ret;
}
