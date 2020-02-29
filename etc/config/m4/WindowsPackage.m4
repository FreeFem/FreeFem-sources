; Creating a FreeFem++ package for Microsoft Windows with Inno Setup
; $Id$
;;  m4 def
;; `DHOSTOS' = HOSTOS
;; `SIZEOFPTR' = SIZEOFPTR
;; ifelse(SIZEOFPTR,64, define(`SUFF64',`-64' ),define(`SUFF64',`' ))
;; define(IFMPI,ifelse(len(MPIPROG),0,; ,))
;; define(IFMGW32,ifelse(SIZEOFPTR,32,,;))
;; define(IFMGW64,ifelse(SIZEOFPTR,64,,;))


;; -- end def m4
; The Inno Setup configuration file WindowsPackage.iss is built from
; WindowsPackage.m4 with the command "make WindowsPackage.iss".

; No source file here. They are in the source tar ball.
; suppress -cs version no fltk to day , wait the next version
;  FH version 3.0-1
[Setup]
AppName=FreeFem++-win`'SIZEOFPTR-VERSION
AppVerName=FreeFem++ version VERSION (win SIZEOFPTR bits)
DefaultDirName={pf}\FreeFem++`'SUFF64
DefaultGroupName=FreeFem++`'SUFF64

Compression=lzma
SolidCompression=yes
ChangesAssociations=yes
OutputBaseFilename=FreeFem++-VERSION-win`'SIZEOFPTR
ChangesEnvironment=yes

[Dirs]
Name: "{app}";
; set writing permissions for  examples with write and  read files
Name: "{app}\examples\misc"; Permissions: everyone-full
Name: "{app}\examples\plugin"; Permissions: everyone-full
Name: "{app}\examples\tutorial"; Permissions: everyone-full
Name: "{app}\examples\3d"; Permissions: everyone-full
Name: "{app}\examples\3dSurf"; Permissions: everyone-full
Name: "{app}\examples\3dCurve"; Permissions: everyone-full
Name: "{app}\examples\examples"; Permissions: everyone-full
Name: "{app}\examples\eigen"; Permissions: everyone-full
Name: "{app}\idp"; Permissions: everyone-full
IFMPI Name: "{app}\examples\mpi"; Permissions: everyone-full
IFMPI Name: "{app}\examples\hpddm"; Permissions: everyone-full


[Files]
; README
Source: "README.md"; DestDir: "{app}"
Source: "readme\README_WINDOWS.md"; DestDir: "{app}"
Source: "readme\INNOVATION"; DestDir: "{app}"
Source: "readme\AUTHORS"; DestDir: "{app}"
Source: "readme\BUGS"; DestDir: "{app}"
Source: "readme\COPYRIGHT"; DestDir: "{app}"
Source: "readme\COPYING"; DestDir: "{app}"
;Source: "README"; DestDir: "{app}"
;Source: "crimson-freefem++.zip"; DestDir: "{app}"
;Source: "0ldUserReadMe.txt"; DestDir: "{app}"

; Programs
Source: "src\bin-win32\FreeFem++.exe"; DestDir: "{app}"
Source: "freefem++.pref"; DestDir: "{app}"
ifelse(len(MPIPROG),0,; ,)Source: "src\bin-win32\FreeFem++-mpi.exe"; DestDir: "{app}"
ifelse(len(MPIPROG),0,; ,)Source: "src\mpi\ff-mpirun"; DestDir: "{app}"
Source: "src\bin-win32\launchff++.exe"; DestDir: "{app}"
;   no freefem++-cs today see ALH (FH)
;Source: "src\bin-win32\FreeFem++-cs.exe"; DestDir: "{app}"
;Source: "src\ide\FreeFem++-cs.exe"; DestDir: "{app}"
Source: "src\nw\ffglut.exe"; DestDir: "{app}"
Source: "src\medit\ffmedit.exe"; DestDir: "{app}"
;Source: "src\bin-win32\FreeFem++-nw.exe"; DestDir: "{app}"
Source: "src\bin-win32\bamg.exe"; DestDir: "{app}"
Source: "src\bin-win32\cvmsh2.exe"; DestDir: "{app}"
; Source: "src\bin-win32\drawbdmesh.exe"; DestDir: "{app}"
Source: "src\bin-win32\*.dll"; DestDir: "{app}"
Source: "plugin\seq\ff-c++"; DestDir: "{app}"
Source: "plugin\seq\ff-get-dep.awk"; DestDir: "{app}"
Source: "plugin\seq\WHERE_LIBRARY-config"; DestDir: "{app}"
Source: "plugin\seq\WHERE_LIBRARY"; DestDir: "{app}"
Source: "plugin\seq\WHERE_LIBRARY-download"; DestDir: "{app}"
Source: "plugin\seq\ff-pkg-download"; DestDir: "{app}"
Source: "plugin\seq\ff-get-dep"; DestDir: "{app}"

; mingwm10.dll is necessary when "-mthreads" is used as a compilation
; flag.
;
;ldd.exe src/bin-win32/FreeFem++.exe  |awk '/mingw64/ {print "cygpath -w ",$3}'|sh|awk '{print "IFMGW64 Source: @" $0 "@ DestDir: @{app}@"}'|sed 's/@/"/g'



IFMGW32 ; mingw32  ....    FH. I have put all dll in bin-win32 dir ....
IFMGW32 Source: "C:\MinGW\bin\mingwm10.dll"; DestDir: "{app}"
; Source: "C:\Cygwin\bin\glut32.dll"; DestDir: "{app}"
IFMGW32 Source: "C:\MinGW\msys\1.0\bin\freeglut.dll"; DestDir: "{app}"
IFMGW32 Source: "C:\MinGW\bin\pthreadGC2.dll"; DestDir: "{app}"
IFMGW32 Source: "C:\MinGW\bin\libgcc_s_dw2-1.dll"; DestDir: "{app}"
IFMGW32 Source: "C:\MinGW\bin\libstdc++-6.dll"; DestDir: "{app}"
IFMGW32 Source: "C:\MinGW\bin\libgfortran-*.dll"; DestDir: "{app}"
IFMGW32 Source: "C:\MinGW\bin\libquadmath-*.dll"; DestDir: "{app}"

IFMGW64 Source: "C:\msys64\mingw64\bin\libgcc_s_seh-*.dll"; DestDir: "{app}"
IFMGW64 Source: "C:\msys64\mingw64\bin\libstdc++-*.dll"; DestDir: "{app}"
IFMGW64 Source: "C:\msys64\mingw64\bin\libwinpthread-1.dll"; DestDir: "{app}"
IFMGW64 Source: "C:\msys64\mingw64\bin\libgfortran-*.dll"; DestDir: "{app}"
IFMGW64 Source: "C:\msys64\mingw64\bin\libquadmath-*.dll"; DestDir: "{app}"
IFMGW64 Source: "C:\msys64\mingw64\bin\libfreeglut.dll"; DestDir: "{app}"


IFMGW64 ; mingw64 ....   FH. I have put all dll in bin-win32 dir ....

;; end of mingw ------------


; Does not include FreeFem++-x11 which would need the Cygwin X-Server
; Does not include FreeFem++-glx which would need the Cygwin X-Server

; Examples
Source: "idp\*.idp"; DestDir: "{app}\idp"
Source: "examples\misc\*.edp"; DestDir: "{app}\examples\misc"
Source: "examples\eigen\*.edp"; DestDir: "{app}\examples\eigen"
Source: "examples\tutorial\*.edp"; DestDir: "{app}\examples\tutorial"
Source: "examples\tutorial\aile.msh"; DestDir: "{app}\examples\tutorial"
Source: "examples\tutorial\xyf"; DestDir: "{app}\examples\tutorial"
Source: "examples\examples\*.edp"; DestDir: "{app}\examples\examples"
Source: "examples\plugin\*.edp"; DestDir: "{app}\examples\plugin"
Source: "examples\plugin\*.pgm"; DestDir: "{app}\examples\plugin"
Source: "examples\plugin\*.pts"; DestDir: "{app}\examples\plugin"
Source: "examples\plugin\cube.msh"; DestDir: "{app}\examples\plugin"
Source: "examples\plugin\g.gmesh"; DestDir: "{app}\examples\plugin"

Source: "plugin\seq\load.link"; DestDir: "{app}\plugin"
Source: "plugin\include-tmp\*"; DestDir: "{app}\include"
Source: "examples\3d\*.edp"; DestDir: "{app}\examples\3d"
Source: "examples\3d\dodecaedre01.mesh"; DestDir: "{app}\examples\3d"
Source: "examples\3d\lac-leman-v4.msh"; DestDir: "{app}\examples\3d"
Source: "examples\3dSurf\*.edp"; DestDir: "{app}\examples\3dSurf"
Source: "examples\3dCurve\*.edp"; DestDir: "{app}\examples\3dCurve"
IFMPI Source: "examples\mpi\ff*.txt"; DestDir: "{app}\examples\mpi"
IFMPI Source: "examples\mpi\*.edp"; DestDir: "{app}\examples\mpi"
IFMPI Source: "examples\hpddm\*.edp"; DestDir: "{app}\examples\hpddm"
;Source: "0ldUserReadMe.txt"; DestDir: "{app}\examples\load"
;Source: "0ldUserReadMe.txt"; DestDir: "{app}\examples\tutorial"
;Source: "0ldUserReadMe.txt"; DestDir: "{app}\examples\examples"
;Source: "0ldUserReadMe.txt"; DestDir: "{app}\examples\misc"
;Source: "0ldUserReadMe.txt"; DestDir: "{app}\examples\eigen"



; Documentation files may need to be copied from another machine if
; Cygwin refuses to build them.

Source: "FreeFEM-documentation.pdf"; DestDir: "{app}"

; Icons for Windows can be created from a 32x32 image with icotool
; (Linux Debian unstable), or IrfanView (Windows, not very good
; results) or paint (Windows, save in .bmp then rename to .ico).

Source: "etc\logo\logo.ico"; DestDir: "{app}"

[Icons]

; Menu
Name: "{group}\FreeFem++"; Filename: "{app}\launchff++.exe"; IconFilename: "{app}\logo.ico"
;Name: "{group}\FreeFem++ GUI"; Filename: "{app}\FreeFem++-cs.exe"
Name: "{group}\PDF manual"; Filename: "{app}\freefem++doc.pdf"
Name: "{group}\Examples\Tutorial"; Filename: "{app}\examples\tutorial"
Name: "{group}\Examples\chapt3"; Filename: "{app}\examples\examples"
Name: "{group}\Examples\load"; Filename: "{app}\examples\plugin"
Name: "{group}\Examples\Main"; Filename: "{app}\examples\misc"
Name: "{group}\Examples\Eigenvalues"; Filename: "{app}\examples\eigen"
Name: "{group}\Examples\3d"; Filename: "{app}\examples\3d"
Name: "{group}\Examples\3dSurf"; Filename: "{app}\examples\3dSurf"
Name: "{group}\Examples\3dCurve"; Filename: "{app}\examples\3dCurve"
IFMPI Name: "{group}\Examples\mpi"; Filename: "{app}\examples\mpi"
IFMPI Name: "{group}\Examples\hpddm"; Filename: "{app}\examples\hpddm"
Name: "{group}\Uninstall FreeFem++ VERSION"; Filename: "{uninstallexe}"

; Desktop
Name: "{userdesktop}\FreeFem++ VERSION"; Filename: "{app}\launchff++.exe"; IconFilename: "{app}\logo.ico"
;Name: "{userdesktop}\FreeFem++ VERSION GUI"; Filename: "{app}\FreeFem++-cs.exe"
Name: "{userdesktop}\FreeFem++ VERSION Examples"; Filename: "{group}\Examples"


[Registry]

; Link .edp file extension to FreeFem++
Root: HKCR; Subkey: ".edp"; ValueType: string; ValueName: ""; ValueData: "FreeFemVERSIONScript"; Flags: uninsdeletevalue
Root: HKCR; Subkey: "FreeFemVERSIONScript"; ValueType: string; ValueName: ""; ValueData: "FreeFem++ Script"; Flags: uninsdeletekey
Root: HKCR; Subkey: "FreeFemVERSIONScript\DefaultIcon"; ValueType: string; ValueName: ""; ValueData: "{app}\logo.ico"
Root: HKCR; Subkey: "FreeFemVERSIONScript\shell\open\command"; ValueType: string; ValueName: ""; ValueData: """{app}\launchff++.exe"" ""%1"""


[Tasks]
Name: modifypath; Description: &Add application directory to your system path (if missing you can have trouble with on-the-fly graphic ) ; Flags:  checkedonce
; unchecked

[Code]
function ModPathDir(): TArrayOfString;
var
			Dir:	TArrayOfString;
		begin
			setArrayLength(Dir, 1)
			Dir[0] := ExpandConstant('{app}');
			Result := Dir;
		end;
#include "bin\modpath.iss"
