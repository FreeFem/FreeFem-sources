<!----------------------------------------------------------------------------------->
<!--- This file is part of FreeFEM.                                               --->
<!--- Laboratoire Jacques-Louis Lions                                             --->
<!--- Sorbonne Université, UMR 7598, Paris, F-75005 France                        --->
<!---                                                                             --->
<!--- FreeFEM is free software: you can redistribute it and/or modify             --->
<!--- it under the terms of the GNU Lesser General Public License as published by --->
<!--- the Free Software Foundation, either version 3 of the License, or           --->
<!--- (at your option) any later version.                                         --->
<!---                                                                             --->
<!--- FreeFEM is distributed in the hope that it will be useful,                  --->
<!--- but WITHOUT ANY WARRANTY; without even the implied warranty of              --->
<!--- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               --->
<!--- GNU Lesser General Public License for more details.                         --->
<!---                                                                             --->
<!--- You should have received a copy of the GNU Lesser General Public License    --->
<!--- along with FreeFEM.  If not, see <http://www.gnu.org/licenses/>.            --->
<!----------------------------------------------------------------------------------->

# How to compile FreeFem++ on Microsoft Windows
_april 2017_

Visit this [web page](http://www.freefem.org/ff++/windows.php)

Bug warning (Windows 64): if you launch FreeFem++ without script by double clicking
on the icon, you get an error (it is due to a bug in GetOpenFileName).

From version 3.52 onwards, the Windows 64 version is built with MPI support and
with the following options
```bash
./configure '--enable-download'
```

Then, execute and follow the instructions.

## To test the MPI usage: in windows terminal (cmd, shell, PowerShell)

Go to the directory `C:\Program Files (86)\FreeFem++\examples\mpi`

To lauch an example:
```bash
mpiexec.exe -np 4 FreeFem++-mpi DDM-Schwarz-Lame-2d.edp -wg
```
or without graphics:
```bash
mpiexec.exe -np 4 FreeFem++-mpi DDM-Schwarz-Lame-2d.edp
```

## To use MPI, you must first install `MSMPI`

Download [MS MPI V7](https://www.microsoft.com/en-us/download/details.aspx?id=49926 ),
and install both `msmpisdk.msi` and `MSMpiSetup.exe`.

Remark:

Under msys2 do not forger to do open `c:\msys64\mingw64.ini` in an editor
and remove `rem` before `set MSYS2_PATH_TYPE=inherit`
