; Creating a FreeFem++ package for Microsoft Windows with Inno Setup
; $Id$

; The Inno Setup configuration file is built from this one with the
; command "make WindowsPackage.iss".

; No source file here. They are in the source tar ball.

[Setup]
AppName=FreeFem++-VERSION
AppVerName=FreeFem++ version VERSION
DefaultDirName={pf}\FreeFem++-VERSION
DefaultGroupName=FreeFem++-VERSION
Compression=lzma
SolidCompression=yes
ChangesAssociations=yes
OutputBaseFilename=FreeFem++-VERSION

[Files]
Source: "src\std\FreeFem++.exe"; DestDir: "{app}"
Source: "src\ide\FreeFem++-cs.exe"; DestDir: "{app}"
Source: "src\ide\FreeFem++-cs-server.exe"; DestDir: "{app}"
Source: "src\nw\FreeFem++-nw.exe"; DestDir: "{app}"
Source: "c:/cygwin/bin/cygwin1.dll"; DestDir: "{app}"
Source: "examples++\*.edp"; DestDir: "{app}\examples++"
Source: "examples++-eigen\*.edp"; DestDir: "{app}\examples++-eigen"
Source: "examples++-tutorial\*.edp"; DestDir: "{app}\examples++-tutorial"
Source: "examples++-other\*.edp"; DestDir: "{app}\examples++-other"
Source: "DOC\manual-full.pdf"; DestDir: "{app}"
Source: "DOC\manual-full.ps"; DestDir: "{app}"
Source: "logo.ico"; DestDir: "{app}"

[Icons]

; Menu
Name: "{group}\FreeFem++"; Filename: "{app}\FreeFem++-cs.exe"; IconFilename: "{app}\logo.ico"
Name: "{group}\PDF manual"; Filename: "{app}\manual-full.pdf"
Name: "{group}\Postscript manual"; Filename: "{app}\manual-full.ps"
Name: "{group}\Examples\Tutorial"; Filename: "{app}\examples++-tutorial"
Name: "{group}\Examples\Main"; Filename: "{app}\examples++"
Name: "{group}\Examples\Eigenvalues"; Filename: "{app}\examples++-eigen"
Name: "{group}\Uninstall FreeFem++ VERSION"; Filename: "{uninstallexe}"

; Desktop
Name: "{userdesktop}\FreeFem++ VERSION"; Filename: "{app}\FreeFem++-cs.exe"; IconFilename: "{app}\logo.ico"

[Registry]

; Link .edp file extension to FreeFem++
Root: HKCR; Subkey: ".edp"; ValueType: string; ValueName: ""; ValueData: "FreeFemVERSIONScript"; Flags: uninsdeletevalue
Root: HKCR; Subkey: "FreeFemVERSIONScript"; ValueType: string; ValueName: ""; ValueData: "FreeFem++ Script"; Flags: uninsdeletekey
Root: HKCR; Subkey: "FreeFemVERSIONScript\DefaultIcon"; ValueType: string; ValueName: ""; ValueData: "{app}\logo.ico"
Root: HKCR; Subkey: "FreeFemVERSIONScript\shell\open\command"; ValueType: string; ValueName: ""; ValueData: """{app}\FreeFem++-cs.exe"" ""%1"""
