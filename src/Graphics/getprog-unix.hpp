#ifndef GETPROG_UNIX_HPP_
#define GETPROG_UNIX_HPP_

#include "mode_open.hpp"
#if _WIN32
#include "ff-win32.cpp"
#else
#include <unistd.h>
#endif

// FFCS: required for chdir() when g++ 4.7 is used
#include <unistd.h>

// FFCS: redirecting visualization output
#include "ffapi.hpp"
#include "strversionnumber.hpp"
#include <string>

extern long mpirank;
extern long verbosity;
extern FILE *ThePlotStream; // Add for new plot. FH oct 2008
// for the environnement variables ...
extern const char *prognamearg;
extern const char *edpfilenamearg;
extern bool waitatend;
extern bool consoleatend;
extern bool echo_edp;
extern bool NoGraphicWindow;

char *Shell_Space(const char *s);

char *Shell_Space(const char *s) {
  const char *c = s;
  int nbspace;
  int i;
  for (i = 0; i < 100000; ++i) {
    if (!s[i]) break;
    else if (isspace(s[i])) ++nbspace;
  }
  if (!(i < 100000)) {
    cerr << " Bug Shell_Space routine " << endl;
    exit(1);
  }

#ifdef _WIN32
  char * p = new char[i + 1 + nbspace];
  char * q = p;
  for (i = 0; i < 100000; ++i) {
    if (!s[i]) break;
    else if (isspace(s[i])) {
      *q ++= '^';
      *q ++= s[i];
    }
    else *q ++= s[i];
  }
#else
  char *p = new char[i + nbspace];
  char *q = p;
  for (i = 0; i < 100000; ++i) {
    if (!s[i]) break;
    else if (isspace(s[i])) {
      *q ++= '\\';
      *q ++= s[i];
    }
    else *q ++= s[i];
  }
#endif
  *q ++= '\0';
  assert(q - p <= i + nbspace);
  return p;
}

extern void (*init_lgparallele)();

// <<getprog>> called by [[file:../lglib/lg.ypp::getprog]]
int getprog(char *fn, int argc, char **argv) {
  waitatend = 0; // init_lgparallele==0; // wait if not parallel
  consoleatend = false; // bug with redirection FH
  int ret = 0;
  *fn = '\0';
#ifdef _WIN32
  const int lsuffix = 4;
#else
  const int lsuffix = 0;
#endif

#ifdef PROG_FFGLUT
  const char *ffglut = PROG_FFGLUT;
#else
  const char *ffglut = "ffglut";
#endif
  const char *progffglut = 0;
  const char *fileglut = 0;
  bool noffglut = true;

  // FFCS - remove the test for noffglut to be able to create pictures
  // in any situation. Even FreeFem++-mpilang needs to send pictures
  // (eg when called in a FreeFem++-server situation by EJS)
  // is the name -nw or -nw.exe -> no graphics

  // if no ffglut do try to launch ffglut by default - april 2017 (FH)
#ifdef PROG_FFGLUT
  noffglut = false;
  NoGraphicWindow = false;
#else
  noffglut = true;
  NoGraphicWindow = true;
#endif
  ffapi::ff_ch2edpdtmpir = false;
  bool ch2edpdir = false;
  if(argc)
    prognamearg = argv[0];

  if (prognamearg) { // FH add to remove ffglut in case of -nw or -nw.exe , mpi ,mpi.exp program FH juin 2014 - april 2017
    int l = strlen(prognamearg);
    if( ((l > 4) && (strcmp("-nw", prognamearg + l - 3) == 0))
      || ((l > 8) && (strcmp("-nw.exe", prognamearg + l - 7) == 0))
      || ((l > 5) && (strcmp("-mpi", prognamearg + l - 4) == 0)) // Correct april 2017 FH
      || ((l > 9) && (strcmp("-mpi.exe", prognamearg + l - 8) == 0)) // Correct april 2017 FH
    ) {
      consoleatend = false;
      noffglut = true;
      NoGraphicWindow = true;
      waitatend = false;
    }
  }
  bool flagnw = 0;
  echo_edp = true;
  ffapi::ff_justcompile = false;
  if(argc) {
    for (int i = 1; i < argc; i++) {
      if (verbosity > 3) cout << " ARGV " << i << " " << argv[i] << " " << fn << endl;
      if (ret == 0 && strcmp(argv[i], "-f") == 0 && i + 1 < argc) {
        strcpy(fn, argv[i + 1]);
        i++;
        edpfilenamearg = argv[i];
        ret = 1;
      }
      else if (strcmp(argv[i], "-v") == 0 && i + 1 < argc) {
        verbosity = atoi(argv[i + 1]);
        i++;
        if (verbosity > 10) cout << " verbosity " << verbosity << endl;
      }
      else if ((strcmp(argv[i], "-nw") == 0) || (strcmp(argv[i], "-ng") == 0)) {// add -ng april 2017
        flagnw = true;
        consoleatend = false;
        noffglut = true;
        NoGraphicWindow = true;
        waitatend = false; // add modif FH. juin 2014 ..
      }
      else if (strcmp(argv[i], "-wg") == 0) {
        noffglut = false;
        NoGraphicWindow = false;
      }

      else if (strcmp(argv[i], "-ne") == 0) // no edp
        echo_edp = false;
      else if (strcmp(argv[i], "-cd") == 0)
        ch2edpdir = true;
      else if (strcmp(argv[i], "-cdtmp") == 0)
        ffapi::ff_ch2edpdtmpir = true;
      else if (strcmp(argv[i], "-jc") == 0) {
        ffapi::ff_justcompile = true;
        waitatend = false;
        consoleatend = false;
        noffglut = true;
        NoGraphicWindow = true;
        waitatend = false;
      }
      else if (strcmp(argv[i], "-ns") == 0) // no script
        echo_edp = false;
      else if (strcmp(argv[i], "-nowait") == 0)
        waitatend = false;
      else if (strcmp(argv[i], "-nc") == 0)
        consoleatend = false;
      else if (strcmp(argv[i], "-log") == 0)
        consoleatend = true;
      else if (strcmp(argv[i], "-wait") == 0)
        waitatend = true;
      else if (strcmp(argv[i], "-fglut") == 0 && i + 1 < argc) {
        fileglut = argv[++i];
        noffglut = true;
      }
      else if (strcmp(argv[i], "-glut") == 0 && i + 1 < argc) {
        progffglut = argv[++i];
        if (flagnw) {
          noffglut = true;
          NoGraphicWindow = true;// if -nw => no graphic in anycase
        }
        else {
        noffglut = false;
        NoGraphicWindow = false;
        }
      }
      else if (strcmp(argv[i], "-gff") == 0 && i + 1 < argc) {
        progffglut = Shell_Space(argv[++i]);
        if (flagnw) { // if -nw => no graphic in anycase
          noffglut = true;
          NoGraphicWindow = true;// if -nw => no graphic in anycase
        }
        else {
          noffglut = false;
          NoGraphicWindow = false;
        }
      }
      else if (strcmp(argv[i], "-?") == 0)
        ret = 2;
      else if (strcmp(argv[i], "-f") == 0 && i + 1 < argc) {
        strcpy(fn, argv[++i]);
        ret = 1;
        edpfilenamearg = argv[i];
        if (verbosity > 4) cout << " fn: " << fn << endl;
      }
      else if (ret == 0) {
        strcpy(fn, argv[i]);
        edpfilenamearg = argv[i];
        ret = 1;
        if (verbosity > 4) cout << " fn: " << fn << endl;
      }
    }
  }
  if (ch2edpdir && edpfilenamearg) {
    int i = 0;
    int l = strlen(edpfilenamearg);
#ifdef _WIN32
    const char sepdir = '\\';
#else
    const char sepdir = '/';
#endif

    for (i = l - 1; i >= 0; i--)
      if (edpfilenamearg[i] == sepdir) break;

    if (i > 0) {
      char *dir = new char[l + 1];
      strcpy(dir, edpfilenamearg);
      dir[i] = 0;
      int err = 0;
      if (verbosity > 1)
        cout << " chdir '" << dir << "'" << endl;
      // FFCS: mingw64 API change
      err = chdir(dir);

      //cout << err << endl;
      if (err) {
        cerr << " error: chdir " << dir << endl;
        exit(1);
      }
      delete[] dir;
    }
  }
  if (!progffglut && !noffglut)
    progffglut = ffglut;

  if (progffglut && mpirank == 0) {
    // FFCS: divert stream to FFCS
    ThePlotStream = ffapi::ff_popen(progffglut, "w");

    // FFCS: just forget the following message because it could be understood as an error during FFCS execution
    // although ffglut is not used there.

    //if(verbosity)
    // printf(" EXEC of the plot: %s\n",progffglut);

    if (!ThePlotStream) { cerr << " Error popen "<< progffglut << endl; exit(1); }
  }
  else if (fileglut && mpirank == 0) { // correction progffglut -> fileglut v3.0-2 FH.
    ThePlotStream = fopen(fileglut, MODE_WRITE_BINARY);
    if (verbosity)
      cout << " save of the plot in file: " << fileglut << endl;
    if (!ThePlotStream) {
      cerr << " Error save file glut " << fileglut
           << " mode " << MODE_WRITE_BINARY << endl;
      exit(1);
    }
  }

#ifdef _WIN32
  if (ret == 0) {
    if (ShowOpenDialogBox1(fn))
      ret = 1;
  }
#endif

  if (ret != 1) {
    const char *ff = argc ? argv[0] : "FreeFem++";

    cout << ff << " - version " << StrVersionNumber() << sizeof(void*)*8 << "bits" << endl;
    cout << "License: LGPL 3+" << endl;
    cout << "Usage: " << ff << " [FreeFEM arguments] filename [script arguments]" << endl;
    cout << "FreeFEM arguments:" << endl;
    cout << "\t-f:     [filename]  script file name" << endl;
    cout << "\t-v:     [verbosity] level of FreeFEM output (0 - 1000000)" << endl;
    cout << "\t-nw:                no graphics" << endl;
    cout << "\t-wg:                with graphics" << endl;
    cout << "\t-ne:                no edp script output" << endl;
    cout << "\t-cd:                change directory to script directory"<< endl;
    cout << "\t-cdtmp:             change directory to /tmp" << endl;
    cout << "\t-jc:                just compile" << endl;
    cout << "\t-ns:                same as -ne" << endl;
    cout << "\t-nowait:            do not wait graphics at the end" << endl;
    cout << "\t-nc:                without console (MS Windows only)" << endl;
    cout << "\t-log:               with console (MS Windows only)" << endl;
    cout << "\t-wait:              wait graphics at the end" << endl;
    cout << "\t-fglut: [filename]  redirect graphics in file" << endl;
    cout << "\t-glut:  [command]   use custom glut" << endl;
    cout << "\t-gff:   [command]   use custom glut (with space quoting)" << endl;
    cout << "\t-?:                 show help" << endl << endl;

    if (noffglut) cout << "without default ffglut: " << ffglut << endl;
    else cout << "with default ffglut: " << ffglut << endl;

    cout << endl;
    cout << "FreeFEM website: https://freefem.org/" << endl;
    cout << "FreeFEM documentation: https://doc.freefem.org/" << endl;
    cout << "FreeFEM forum: https://community.freefem.org/" << endl;
    cout << "FreeFEM modules: https://modules.freefem.org/" << endl;
    exit(1);
    return ret;
  }
  if (verbosity > 10)
    cout << " file: " << fn << endl;
  // cout << " verbosity="<< verbosity << endl;
  return 1;
}

#endif //GETPROG_UNIX_HPP_
