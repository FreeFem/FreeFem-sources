/****************************************************************************/
/* This file is part of FreeFEM.                                            */
/*                                                                          */
/* FreeFEM is free software: you can redistribute it and/or modify          */
/* it under the terms of the GNU Lesser General Public License as           */
/* published by the Free Software Foundation, either version 3 of           */
/* the License, or (at your option) any later version.                      */
/*                                                                          */
/* FreeFEM is distributed in the hope that it will be useful,               */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of           */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            */
/* GNU Lesser General Public License for more details.                      */
/*                                                                          */
/* You should have received a copy of the GNU Lesser General Public License */
/* along with FreeFEM. If not, see <http://www.gnu.org/licenses/>.          */
/****************************************************************************/
// SUMMARY : ...
// LICENSE : LGPLv3
// ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
// AUTHORS : Frederic Hecht
// E-MAIL  : frederic.hecht@sorbonne-universite.fr

#define NOMINMAX 1
#define byte win_byte_override
#include <windows.h>

#include <direct.h>
#include <cstring>
#include <cstdlib>
#include <string>
#include <cstring>
#include <iostream>
#include <cassert>
using namespace std;

const char C = '"';

/*!
 *
 */
char *dirname (const char *pp) {
  int i, l = strlen(pp);
  for (i = l-1; i >= 0; i--)
    if (pp[i] == '\\') break;
  char *dir = new char[l+1];
  strcpy(dir, pp);
  dir[i] = 0;
  return dir;
}

/*!
 *
 */
char *newq (const char *p) {
  int l = strlen(p);
  assert(l < 1000);
  char *c = new char[l+100];
  assert(c);
  c[0] = C;
  strcpy(c+1, p);
  c[l+1] = C;
  c[l+2] = 0;
  return c;
}

/*!
 *
 */
int main (int argc, const char **argv) {
  char filename[MAX_PATH];
  const int n100 = 100;
  const char *ffargv[n100];
  int debug = 0;
  char *dir = 0;
  const char *dirll = dirname(argv[0]);
  const char *pp = 0;
  char ff[1000];
  const char *fe = "\\freefem++.exe";//"C:\\msys64\\home\\hecht\\ff++\\src\\bin-win32\\ar v.exe";
  *ff = 0;
  strcat(ff, dirll);
  strcat(ff, fe);
  //string cmd="freefem++.exe ";
  for (int i = 0; i < n100; ++i)
    ffargv[i] = 0;
  int nffa = 0;
  ffargv[nffa++] = newq(ff);

  for (int i = 1; i < argc; ++i) {
    if (strcmp("++d", argv[i]) == 0)
      debug = 1;
    else {
      // cmd += C;
      // cmd += argv[i];
      if (!pp && strlen(argv[i]) > 2)
        if (argv[i][1] == ':')
          pp = argv[i];
      // cmd += C;
      // cmd += " ";
      ffargv[nffa++] = newq(argv[i]);

      if( debug) cout << "  ffl: arg " << i << argv[i] << endl;
    }
  }
  if(!pp) {
    //      cerr << " Sorry no file name "<< endl;
    //      cerr << " Drag and Drop the file icon on the application icon or double clip on script file" << endl;
    //   cmd += " -wait -log";

    OPENFILENAME ofn;
    ZeroMemory(&filename, sizeof(filename));
    ZeroMemory(&ofn, sizeof(ofn));
    ofn.lStructSize = sizeof( ofn );
    ofn.hwndOwner   = NULL; // If you have a window to center over, put its HANDLE here
    ofn.lpstrFilter = "freefem++ Files (*.edp)\0*.edp\0All Files (*.*)\0*.*\0\0";
    ofn.lpstrFile   = filename;
    ofn.nMaxFile    = MAX_PATH;
    ofn.lpstrTitle  = "Please, select a file";
    ofn.Flags       = OFN_DONTADDTORECENT | OFN_FILEMUSTEXIST;

    if (GetOpenFileNameA(&ofn)) {
      std::cout << "You chose the file \"" << filename << "\"\n";
      /*	cmd +=" -cd ";
      //cmd += C;
      for(int i=0; i< MAX_PATH; ++i)
      if( filename[i] == 0 ) break;
      else if( filename[i] == ' ' ) cmd += "\\ ";
      else cmd += filename[i] ;
      */
      if (!pp && strlen(filename) > 2)
        if (filename[1] == ':')
          pp = filename;
      ffargv[nffa++] = newq(filename);
      //cmd += filename;
      //cmd += C;
      // cout << " system : " << cmd << endl;
      //int ret= system(cmd.c_str());
    } else {
      // All this stuff below is to tell you exactly how you messed up above.
      // Once you've got that fixed, you can often (not always!) reduce it to a 'user cancelled' assumption.
      switch (CommDlgExtendedError()) {
        case CDERR_DIALOGFAILURE   : std::cout << "CDERR_DIALOGFAILURE\n";   break;
        case CDERR_FINDRESFAILURE  : std::cout << "CDERR_FINDRESFAILURE\n";  break;
        case CDERR_INITIALIZATION  : std::cout << "CDERR_INITIALIZATION\n";  break;
        case CDERR_LOADRESFAILURE  : std::cout << "CDERR_LOADRESFAILURE\n";  break;
        case CDERR_LOADSTRFAILURE  : std::cout << "CDERR_LOADSTRFAILURE\n";  break;
        case CDERR_LOCKRESFAILURE  : std::cout << "CDERR_LOCKRESFAILURE\n";  break;
        case CDERR_MEMALLOCFAILURE : std::cout << "CDERR_MEMALLOCFAILURE\n"; break;
        case CDERR_MEMLOCKFAILURE  : std::cout << "CDERR_MEMLOCKFAILURE\n";  break;
        case CDERR_NOHINSTANCE     : std::cout << "CDERR_NOHINSTANCE\n";     break;
        case CDERR_NOHOOK          : std::cout << "CDERR_NOHOOK\n";          break;
        case CDERR_NOTEMPLATE      : std::cout << "CDERR_NOTEMPLATE\n";      break;
        case CDERR_STRUCTSIZE      : std::cout << "CDERR_STRUCTSIZE\n";      break;
        case FNERR_BUFFERTOOSMALL  : std::cout << "FNERR_BUFFERTOOSMALL\n";  break;
        case FNERR_INVALIDFILENAME : std::cout << "FNERR_INVALIDFILENAME\n"; break;
        case FNERR_SUBCLASSFAILURE : std::cout << "FNERR_SUBCLASSFAILURE\n"; break;
        default                    : std::cout << "You cancelled.\n";
      }
      //int ret= system(cmd.c_str());
      return 1;
    }
  }
  debug = 0;
  if (pp) {
    if (debug) cout << "  ffl: file:" << pp << endl;
    int i = 0;
    dir = dirname(pp);
    if (debug)
      cout << "  ffl: chdir to " << dir << endl;
    _chdir(dir);
    delete[] dir;
  }
  ffargv[nffa++] = "-log";
  ffargv[nffa++] = "-wait";

  //  cmd += " -wait -log";
  if (debug) {
    // cout << "exec " << cmd << endl;
    for (int i = 0; i < n100; ++i)
      if (ffargv[i])
        cout << i << " " << ffargv[i] << endl;
      else break;
    cout << " call "<<endl;
  }

  int ret = _execvp(ff, ffargv);
  // int ret= system(cmd.c_str());
  return ret;
}
