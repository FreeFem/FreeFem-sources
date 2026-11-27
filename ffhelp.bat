@echo off
rem ffhelp: the online help of FreeFEM from the Windows command line (cmd, PowerShell), as the shell script ffhelp
rem   ffhelp                     the list of the topics with a written help
rem   ffhelp name [name ...]     the help of a function, a type or a keyword
rem   ffhelp -l plugin name      load a plugin first (ex. ffhelp -l medit medit); -l can be repeated
rem   ffhelp -mpi name           with FreeFem++-mpi.exe
rem the FreeFem++.exe of the directory of this file is used (else the one of the PATH)
setlocal EnableExtensions
set "FFDIR=%~dp0"
rem (the directory of this file: before shift, which also shifts %0)
set "FFEXE=FreeFem++.exe"
set "TMPD=%TEMP%\ffhelp_%RANDOM%%RANDOM%"
mkdir "%TMPD%" || exit /b 1
set "LOADS=%TMPD%\loads.txt"
set "NAMES=%TMPD%\names.txt"
type nul > "%LOADS%"
type nul > "%NAMES%"

:args
if "%~1"=="" goto run
if /i "%~1"=="-h" goto usage
if /i "%~1"=="--help" goto usage
if /i "%~1"=="-mpi" goto optmpi
if /i "%~1"=="-l" goto optload
>>"%NAMES%" echo help("%~1");
shift
goto args
:optmpi
set "FFEXE=FreeFem++-mpi.exe"
shift
goto args
:optload
if "%~2"=="" (echo ffhelp: -l needs a plugin name 1>&2 & goto fail)
>>"%LOADS%" echo load "%~2"
shift
shift
goto args

:run
set "FF=%FFDIR%%FFEXE%"
if not exist "%FF%" set "FF=%FFEXE%"
for %%A in ("%NAMES%") do if %%~zA==0 >>"%NAMES%" echo help();
copy /b "%LOADS%"+"%NAMES%" "%TMPD%\ffhelp.edp" > nul
pushd "%TMPD%"
"%FF%" -nw -v 0 ffhelp.edp
set "RC=%ERRORLEVEL%"
popd
rmdir /s /q "%TMPD%"
exit /b %RC%

:usage
echo ffhelp: the online help of FreeFEM from the Windows command line (cmd, PowerShell)
echo   ffhelp                     the list of the topics with a written help
echo   ffhelp name [name ...]     the help of a function, a type or a keyword
echo   ffhelp -l plugin name      load a plugin first (ex. ffhelp -l medit medit); -l can be repeated
echo   ffhelp -mpi name           with FreeFem++-mpi.exe
echo the FreeFem++.exe of the directory of this file is used (else the one of the PATH)
rmdir /s /q "%TMPD%"
exit /b 0
:fail
rmdir /s /q "%TMPD%"
exit /b 1
