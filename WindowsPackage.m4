; Creating a FreeFem++ package for Microsoft Windows with Inno Setup
; $Id$

; The Inno Setup configuration file WindowsPackage.iss is built from
; WindowsPackage.m4 with the command "make WindowsPackage.iss".

; No source file here. They are in the source tar ball.
; suppress -cs version no fltk to day , wait the next version
;  FH version 3.0-1
[Setup]
AppName=FreeFem++-VERSION
AppVerName=FreeFem++ version VERSION
DefaultDirName={pf}\FreeFem++
DefaultGroupName=FreeFem++
Compression=lzma
SolidCompression=yes
ChangesAssociations=yes
OutputBaseFilename=FreeFem++-VERSION
ChangesEnvironment=yes

[Files]
; README 
Source: "README"; DestDir: "{app}"
Source: "README_WINDOWS"; DestDir: "{app}"
Source: "INNOVATION"; DestDir: "{app}"
Source: "AUTHORS"; DestDir: "{app}"
Source: "BUGS"; DestDir: "{app}"
Source: "COPYRIGHT"; DestDir: "{app}"
Source: "COPYING"; DestDir: "{app}"
Source: "README"; DestDir: "{app}"
Source: "crimson-freefem++.zip"; DestDir: "{app}"
Source: "0ldUserReadMe.txt"; DestDir: "{app}"

; Programs
Source: "src\bin-win32\FreeFem++.exe"; DestDir: "{app}"
Source: "src\bin-win32\launchff++.exe"; DestDir: "{app}"
;  to day the dll version do not works so we use the static one (FH)
;Source: "src\bin-win32\FreeFem++-cs.exe"; DestDir: "{app}"
;Source: "src\ide\FreeFem++-cs.exe"; DestDir: "{app}"
Source: "src\nw\ffglut.exe"; DestDir: "{app}"
Source: "src\medit\ffmedit.exe"; DestDir: "{app}"
Source: "src\bin-win32\FreeFem++-nw.exe"; DestDir: "{app}"
Source: "src\bin-win32\bamg.exe"; DestDir: "{app}"
Source: "src\bin-win32\cvmsh2.exe"; DestDir: "{app}"
Source: "src\bin-win32\drawbdmesh.exe"; DestDir: "{app}"
Source: "src\bin-win32\*.dll"; DestDir: "{app}"
Source: "examples++-load\ff-c++"; DestDir: "{app}"

; mingwm10.dll is necessary when "-mthreads" is used as a compilation
; flag.

Source: "C:\Cygwin\bin\mingwm10.dll"; DestDir: "{app}"
Source: "C:\Cygwin\bin\glut32.dll"; DestDir: "{app}"

; Does not include FreeFem++-x11 which would need the Cygwin X-Server
; Does not include FreeFem++-glx which would need the Cygwin X-Server

; Examples
Source: "examples++\*.edp"; DestDir: "{app}\examples++"
Source: "examples++-eigen\*.edp"; DestDir: "{app}\examples++-eigen"
Source: "examples++-tutorial\*.edp"; DestDir: "{app}\examples++-tutorial"
Source: "examples++-tutorial\*.idp"; DestDir: "{app}\examples++-tutorial"
Source: "examples++-tutorial\aile.msh"; DestDir: "{app}\examples++-tutorial"
Source: "examples++-tutorial\xyf"; DestDir: "{app}\examples++-tutorial"
Source: "examples++-chapt3\*.edp"; DestDir: "{app}\examples++-chapt3"
Source: "examples++-other\*.edp"; DestDir: "{app}\examples++-other"
Source: "examples++-load\*.edp"; DestDir: "{app}\examples++-load"
Source: "examples++-load\*.cpp"; DestDir: "{app}\examples++-load"
Source: "examples++-load\load.link"; DestDir: "{app}\examples++-load"
Source: "examples++-load\include-tmp\*"; DestDir: "{app}\examples++-load\include"
Source: "examples++-3d\*.edp"; DestDir: "{app}\examples++-3d"
Source: "examples++-3d\dodecaedre01.mesh"; DestDir: "{app}\examples++-3d"
Source: "examples++-3d\lac-leman-v4.msh"; DestDir: "{app}\examples++-3d"
Source: "examples++-3d\*.idp"; DestDir: "{app}\examples++-3d"
Source: "0ldUserReadMe.txt"; DestDir: "{app}\examples++-load"
Source: "0ldUserReadMe.txt"; DestDir: "{app}\examples++-tutorial"
Source: "0ldUserReadMe.txt"; DestDir: "{app}\examples++-chapt3"
Source: "0ldUserReadMe.txt"; DestDir: "{app}\examples++"
Source: "0ldUserReadMe.txt"; DestDir: "{app}\examples++-eigen"



; Documentation files may need to be copied from another machine if
; Cygwin refuses to build them.

Source: "DOC\freefem++doc.pdf"; DestDir: "{app}"

; Icons for Windows can be created from a 32x32 image with icotool
; (Linux Debian unstable), or IrfanView (Windows, not very good
; results) or paint (Windows, save in .bmp then rename to .ico).

Source: "logo.ico"; DestDir: "{app}"

[Icons]

; Menu
Name: "{group}\FreeFem++"; Filename: "{app}\launchff++.exe"; IconFilename: "{app}\logo.ico"
;Name: "{group}\FreeFem++ GUI"; Filename: "{app}\FreeFem++-cs.exe"
Name: "{group}\PDF manual"; Filename: "{app}\freefem++doc.pdf"
Name: "{group}\Examples\Tutorial"; Filename: "{app}\examples++-tutorial"
Name: "{group}\Examples\chapt3"; Filename: "{app}\examples++-chapt3"
Name: "{group}\Examples\load"; Filename: "{app}\examples++-load"
Name: "{group}\Examples\Main"; Filename: "{app}\examples++"
Name: "{group}\Examples\Eigenvalues"; Filename: "{app}\examples++-eigen"
Name: "{group}\Examples\3d"; Filename: "{app}\examples++-3d"
Name: "{group}\Uninstall FreeFem++ VERSION"; Filename: "{uninstallexe}"

; Desktop
Name: "{userdesktop}\FreeFem++ VERSION"; Filename: "{app}\FreeFem++.exe"; IconFilename: "{app}\logo.ico"
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
 #include "modpath.iss"

