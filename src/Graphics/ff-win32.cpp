
#ifdef WIN32
#include <string>
using namespace std;
#include <windows.h>
#include <commdlg.h>
#include <io.h>      //*OT  use for the console window

BOOL ShowOpenDialogBox1(char *fileName)
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



string ChangeExt(string & ff,const char * suff)
{
  int dot = ff.rfind(".edp");
  assert(dot>0);
  return ff.substr(0,dot)+suff;
}

bool GetConsoleBuff(const string &edpname)
{
  CONSOLE_SCREEN_BUFFER_INFO csbi; //* to get buffer info
  HANDLE hConOut= GetStdHandle(STD_OUTPUT_HANDLE);
  GetConsoleScreenBufferInfo(hConOut, &csbi);
  
  COORD coordLine = {0,0};
  CHAR *szLine=0;  //* buffer to read from the console (a line)
  DWORD dwCharsRead;
  FILE *fp;
  string  fname=ChangeExt(edpname,"log");
  if ((fp = fopen(fname.c_str(),"w"))==NULL) {
    perror(fname.c_str());
    return false;
  }
  
  szLine = new CHAR [csbi.dwSize.X+1];
  for (int i=0; i<csbi.dwCursorPosition.Y; i++) 
	{
        if (ReadConsoleOutputCharacter(hConOut, szLine,
				   csbi.dwSize.X, coordLine,
				   &dwCharsRead)== FALSE)
      {
	perror("ReadConsoleOutputCharacter");
	return false;
      }
    int j=csbi.dwSize.X-1;
    while ((szLine[j] == ' ') && (j > 0)) szLine[j--] =0;
    if (j < csbi.dwSize.X-1) szLine[j+1] = '\n'; 
    fprintf(fp,"%s",szLine);
    coordLine.Y++;
  }
  fclose(fp);
  delete [] szLine;
  return true;
}
#endif
