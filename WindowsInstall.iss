; Installation of FreeFem++ on Microsoft Windows with Inno Setup
; $Id$

[Setup]
AppName=FreeFem++
AppVerName=FreeFem++ version 1.41
DefaultDirName={pf}\FreeFem++
DefaultGroupName=FreeFem++
Compression=lzma
SolidCompression=yes
ChangesAssociations=yes

[Files]
Source: "*"; DestDir: "{app}"; Excludes: "*.o,*.a,*.Po,CVS"; Flags: recursesubdirs

[Icons]
Name: "{group}\FreeFem++"; Filename: "{app}\src\std\FreeFem++.exe"; IconFilename: "{app}\logo.ico"
Name: "{userdesktop}\FreeFem++"; Filename: "{app}\src\std\FreeFem++.exe"; IconFilename: "{app}\logo.ico"
Name: "{group}\examples"; Filename: "{app}\examples++"
Name: "{group}\Uninstall FreeFem++"; Filename: "{uninstallexe}"

[Registry]

; Link .edp file extension to FreeFem++
Root: HKCR; Subkey: ".edp"; ValueType: string; ValueName: ""; ValueData: "FreeFemScript"; Flags: uninsdeletevalue
Root: HKCR; Subkey: "FreeFemScript"; ValueType: string; ValueName: ""; ValueData: "FreeFem++ Script"; Flags: uninsdeletekey
Root: HKCR; Subkey: "FreeFemScript\DefaultIcon"; ValueType: string; ValueName: ""; ValueData: "{app}\logo.ico"
Root: HKCR; Subkey: "FreeFemScript\shell\open\command"; ValueType: string; ValueName: ""; ValueData: """{app}\src\std\FreeFem++.exe"" ""%1"""

