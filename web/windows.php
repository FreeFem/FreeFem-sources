<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01//EN" "http://www.w3.org/TR/html4/strict.dtd">
<html>
<head>
  <meta http-equiv="Content-Type" content="text/html; charset=utf-8">
  <meta http-equiv="Content-Style-Type" content="text/css">
  <title></title>
  <meta name="Generator" content="Cocoa HTML Writer">
  <meta name="CocoaVersion" content="1138.51">
  <style type="text/css">
    p.p1 {margin: 0.0px 0.0px 0.0px 0.0px; font: 14.0px Helvetica}
    p.p2 {margin: 0.0px 0.0px 0.0px 0.0px; font: 19.0px Helvetica; min-height: 23.0px}
    p.p3 {margin: 0.0px 0.0px 0.0px 0.0px; font: 15.0px Helvetica}
    p.p4 {margin: 0.0px 0.0px 0.0px 0.0px; font: 15.0px Helvetica; min-height: 18.0px}
    p.p5 {margin: 0.0px 0.0px 0.0px 0.0px; font: 15.0px Helvetica; color: #0000ee}
    p.p6 {margin: 0.0px 0.0px 0.0px 0.0px; font: 12.0px Helvetica}
    p.p7 {margin: 0.0px 0.0px 0.0px 0.0px; font: 15.0px Courier; color: #008e53}
    p.p8 {margin: 0.0px 0.0px 0.0px 0.0px; font: 16.0px Courier; color: #008e53}
    p.p9 {margin: 0.0px 0.0px 0.0px 0.0px; font: 16.0px Arial}
    p.p10 {margin: 0.0px 0.0px 0.0px 0.0px; font: 15.0px Courier; min-height: 18.0px}
    p.p11 {margin: 0.0px 0.0px 0.0px 0.0px; font: 15.0px Helvetica; color: #ff0000}
    p.p12 {margin: 0.0px 0.0px 0.0px 0.0px; font: 15.0px Helvetica; color: #3e00ff}
    p.p13 {margin: 0.0px 0.0px 0.0px 0.0px; font: 15.0px Helvetica; color: #3e00ff; min-height: 18.0px}
    p.p14 {margin: 0.0px 0.0px 0.0px 0.0px; font: 12.0px Helvetica; min-height: 14.0px}
    p.p15 {margin: 0.0px 0.0px 0.0px 0.0px; font: 15.0px Helvetica; color: #3e1bfe}
    p.p16 {margin: 0.0px 0.0px 0.0px 0.0px; font: 15.0px Helvetica; color: #3e02ff}
    p.p17 {margin: 0.0px 0.0px 0.0px 0.0px; font: 15.0px Helvetica; color: #a20092}
    p.p18 {margin: 0.0px 0.0px 0.0px 0.0px; font: 15.0px Helvetica; color: #1c42d9}
    p.p19 {margin: 0.0px 0.0px 0.0px 0.0px; font: 15.0px Helvetica; color: #6c007d; min-height: 18.0px}
    p.p20 {margin: 0.0px 0.0px 0.0px 0.0px; font: 15.0px Helvetica; color: #6c007d}
    span.s1 {font: 19.0px Helvetica}
    span.s2 {color: #ff2a1a}
    span.s3 {color: #008e53}
    span.s4 {color: #000000}
    span.s5 {text-decoration: underline}
    span.s6 {color: #3e02ff}
    span.s7 {font: 15.0px Helvetica; color: #000000}
  </style>
</head>
<body>
<?php
//  comment in phpp ..
include 'phpfiles.php';
?>
f
<p class="p1">To install the precompile windows package just<span class="Apple-converted-space"> </span></p>
<p class="p1">download the last version from </p>
<p class="p1">Take the file form and download </p> 

    <? pfile($wfile,"Windows 32bits") ?>
    <? pfile($w64file,"Windows 64bits (in test)") ?>
 <p class="p1"> And  execute<span class="Apple-converted-space">  </span>and<span class="Apple-converted-space">  </span>follow the instruction<span class="s1">.</span></p>
<p class="p2"><br></p>
<p class="p2"><span class="Apple-converted-space"> </span></p>
<p class="p3"><span class="s1">-- </span>How to compile FreeFem++ on Microsoft Windows (win32)</p>
<p class="p1"><span class="s1"><span class="Apple-converted-space">    </span></span>F. Hecht<span class="Apple-converted-space">  </span>(Paris, Sept. the 4th, 2013)<span class="Apple-converted-space"> </span></p>
<p class="p3">---------------------------------------------</p>
<p class="p3"><span class="s2">WARNING<span class="Apple-converted-space">  </span>NOW the window version<span class="Apple-converted-space">  </span>is compiled under MINGW </span>for version<span class="Apple-converted-space">  </span>before version 3.20 , oct 7h 2012 , see the end of the file (obsolete now)</p>
<p class="p4"><br></p>
<p class="p3">typo remark: all line in <span class="s3">green</span><span class="Apple-converted-space">  </span>are shell command under mingw32 shell</p>
<p class="p4"><br></p>
<p class="p3">1. Download and install MINGW32<span class="Apple-converted-space"> </span></p>
<p class="p5"><span class="s4">form <a href="http://sourceforge.net/projects/mingw/files/Installer/mingw-get-setup.exe/download">
<span class="s5">http://sourceforge.net/projects/mingw/files/Installer/mingw-get-setup.exe/download</span></a></span></p>
<p class="p6">Answer question in following  windows : </p>
<p class="p6">  7) select  components </p>
<p class="p6">     all basic setup </p>
<p class="p6">     all package : in mingw32  mingw-devlopper-tool, autoconf , automake,  compiler, wget , ... </p>

<p class="p6">  8)   Installation -> Apply changes  </p>
<p class="p4"><br></p>
<p class="p3">2. Under mingw32 shell install wget and unzip</p>
<p class="p7"><span class="Apple-converted-space"> </span>mingw-get install msys-wget</p>
<p class="p7"><span class="Apple-converted-space"> </span>mingw-get.exe install msys-unzip</p>
<p class="p4"><br></p>
<p class="p3">3. To install freeglut of win32 for the graphics part</p>
<p class="p4"><br></p>
<p class="p7">wget http://files.transmissionzero.co.uk/software/development/GLUT/freeglut-MinGW.zip</>
<p class="p7">unzip freeglut-MinGW-2.8.0-1.mp.zip</p>
<p class="p7">cp freeglut/include/GL/* /c/MinGW/include/GL/.</p>
<p class="p7">cp freeglut/lib/lib*.a /c/MinGW/lib/.</p>
<p class="p7">cp freeglut/bin/freeglut.dll /c/MinGW/bin</p>
<p class="p4"><br></p>
<p class="p3">4. install a good blas (OpenBlas) http://xianyi.github.com/OpenBLAS/</p>
<p class="p5"><span class="s4"><span class="Apple-converted-space">   </span>get from<span class="Apple-converted-space">  </span><a href="http://github.com/xianyi/OpenBLAS/tarball/v0.2.3"><span class="s5">http://github.com/xianyi/OpenBLAS/tarball/v0.2.3</span></a></span></p>
<p></p><p></p>
<p class="p7">wget  http://github.com/xianyi/OpenBLAS/tarball/v0.2.3 -O OpenBlas.tgz </p>
<p class="p7">tar zxvf   OpenBlas.tgz </p>
<p class="p7">cd xianyi-OpenBLAS-*/.   </p>
<p class="p7">make </p>
<p class="p7">make install PREFIX=$HOME/soft</p>
<p class="p7">mkdir $HOME/soft/bin
<p class="p7">cp *.dll $HOME/soft/bin
<p></p><p></p>
<p class="p3">5. install MPI for // HPC Pack 2008 R2 Client Pack 4 </p>
<p class="p3">   and  install MPI for // HPC Pack 2008 R2  Pack 4 </p>

<p class="p4"><span class="Apple-converted-space">  </span></p>
<p class="p3">6. install inno setup to build installer :<span class="Apple-converted-space"> </span></p>
<p class="p5"><span class="s4"><span class="Apple-converted-space">  </span><a href="http://www.xs4all.nl/~mlaan2/ispack/isetup-5.4.0.exe"><span class="s5">http://www.xs4all.nl/~mlaan2/ispack/isetup-5.4.0.exe</span></a></span></p>
<p class="p4"><br></p>
<p class="p3">7. GSL for gsl interface is take form</p>
<p class="p5"><span class="s4"><span class="Apple-converted-space">   </span><a href="http://sourceforge.net/projects/mingw-cross/files/%5BLIB%5D%20GSL/mingw32-gsl-1.14-1/"><span class="s5">http://sourceforge.net/projects/mingw-cross/files/%5BLIB%5D%20GSL/mingw32-gsl-1.14-1/</span></a></span></p>
<p class="p4"><br></p>
<p class="p3">9) To download the latest freefem++   tar.gz file contening source form  </p>
    <? pfile($sfile,"Source code") ?>
 <li>  or you can get the latest source from
an anonymous <A href="http://mercurial.selenic.com/" > Mercurial SCM </A> copy with the following unix shell commands&nbsp;:
<p class="p8">hg clone<span class="Apple-converted-space">  </span>http://www.freefem.org/ff++/ff++</p>
<p class="p9">to update do to the last version:</p>
<p class="p8">hg pull</p>
<p class="p8">hg up<span class="Apple-converted-space"> </span></p>
<p class="p4"><br></p>
<p class="p3">To restore, all files build by autoreconf -i command  (automake):</b></p>
<p class="p7">tar zxvf AutoGeneratedFile.tar.gz </p>

<p class="p3">Finaly, the configure argument are:</b></p>
<p class="p3">10)<span class="Apple-converted-space">  </span>Finaly, the configure argument are :</p>
<p class="p7">./configure ’--enable-download’ ’FC=mingw32-gfortran’ ’F77=mingw32-gfortran’ ’CC=mingw32-gcc’ ’CXX=mingw32-g++’ ’-with-blas=$HOME/soft/bin/libopenblas.dll’ ’CXXFLAGS=-I$HOME/soft/include’ ’--enable-generic’ ’--with-wget=wget’ ’MPIRUN=/c/Program Files/Microsoft HPC Pack 2008 R2/Bin/mpiexec.exe</p>
<p class="p10"><br></p>
<p class="p3">-----------------------------------------------------------------------------------------------------------------------------</p>
<p class="p3">Ok until version 3.19-1 (but now this soft is to old to get form the web).<span class="Apple-converted-space"> </span></p>
<p class="p3">FIle version 30/11/2011 F. Hecht.<span class="Apple-converted-space"> </span></p>
<p class="p4"><br></p>
<p class="p11">WARNING<span class="Apple-converted-space">  </span>NOW the window version<span class="Apple-converted-space">  </span>is compiled under MINGW<span class="Apple-converted-space">  </span>(from version 3.11<span class="Apple-converted-space">  </span>14/01/2011 FH)</p>
<p class="p3">So the old dll are incompatible with the new version.<span class="Apple-converted-space"> </span></p>
<p class="p3">It is the fortran compiler under cygwin which is too old<span class="Apple-converted-space">  </span>(not f90 under cygwin).<span class="Apple-converted-space"> </span></p>
<p class="p3">----------------------------------------------------------</p>
<p class="p3">The<span class="Apple-converted-space">  </span>tools to<span class="Apple-converted-space">  </span>be installed are:</p>
<p class="p4"><br></p>
<p class="p3"><b>1) Download and install MINGW32<span class="Apple-converted-space"> </span></b></p>
<p class="p4"><br></p>
<p class="p12">http://sunet.dl.sourceforge.net/project/mingw/Automated%20MinGW%20Installer/mingw-get-inst/mingw-get-inst-20101030/mingw-get-inst-20101030.exe</p>
<p class="p13"><br></p>
<p class="p6">launch:<span class="Apple-converted-space"> </span></p>
<p class="p6"><span class="Apple-converted-space">   </span> mingw-get-inst-20101030.exe</p>
<p class="p14"><br></p>
<p class="p6">Answer question in following  windows : </p>
<p class="p6">   1) do next </p>
<p class="p6">   2) do next </p>
<p class="p6">   3) use preload case </p>
<p class="p6">   4)  accept </p>
<p class="p6">   5 ) select   location of mingw on disk. </p>
<p class="p6">   6) mingw menu name </p>
<p class="p6">  7) select  components</p>
<p class="p6">     all except  ada </p>
<p class="p6">  8)   do install </p>
<p class="p13"><br></p>
<p class="p4"><br></p>
<p class="p3"><b>2) Download and install wget</b> for --enable-download in configure</p>
<p class="p12">http://puzzle.dl.sourceforge.net/project/mingw/mingwPORT/Current%20Releases/wget-1.9.1-mingwPORT.tar.bz2</p>
<p class="p3">under mingw32 shell<span class="Apple-converted-space"> </span></p>
<p class="p3">do:<span class="Apple-converted-space"> </span></p>
<p class="p3">cd /c/users/loginname/download/</p>
<p class="p15"><span class="s4">tar jxvf </span>wget-1.9.1-mingwPORT.tar.bz2</p>
<p class="p15">cp<span class="Apple-converted-space"> </span></p>
<p class="p3">extract and move wget.exe in /usr/bin<span class="Apple-converted-space"> </span></p>
<p class="p4"><br></p>
<p class="p3"><b>3) The glut of win32 from</b></p>
<p class="p4"><br></p>
<p class="p12"><span class="s6">wget<span class="Apple-converted-space">  </span></span>http://web.cs.wpi.edu/~gogo/courses/mingw/winglut.zip</p>
<p class="p12">or<span class="Apple-converted-space"> </span></p>
<p class="p16">wget http://files.transmissionzero.co.uk/software/development/GLUT/freeglut-MinGW.zip</p>
<p class="p4"><br></p>
<p class="p3">The location of include file must be<span class="Apple-converted-space"> </span></p>
<p class="p4"><br></p>
<p class="p17">c:\mingw\include\GL\glut.h</p>
<p class="p17">c:\mingw\include\GL\gl.h</p>
<p class="p17">c:\mingw\include|GL/glu.h</p>
<p class="p4"><br></p>
<p class="p3">add the glut32.dll or freeglut.dll in you directory in the 2 directories:</p>
<p class="p4"><br></p>
<p class="p3">$ find /c/MinGW -name glut</p>
<p class="p17">/c/MinGW/bin/glut32.dll</p>
<p class="p17">/c/MinGW/lib/glut32.dll</p>
<p class="p4"><span class="Apple-converted-space"> </span></p>
<p class="p3"><b>4) the good blas now is:</b><span class="Apple-converted-space"> </span></p>
<p class="p4"><br></p>
<p class="p12">http://www.tacc.utexas.edu/tacc-projects/gotoblas2/downloads/</p>
<p class="p4"><br></p>
<p class="p3">Try to compile<span class="Apple-converted-space"> </span></p>
<p class="p4"><br></p>
<p class="p3"><b>5) install MPI for // versio</b>n<span class="Apple-converted-space"> </span></p>
<p class="p3">HPC Pack 2008 SDK</p>
<p class="p18">http://www.microsoft.com/download/en/details.aspx?id=10505</p>
<p class="p3">HPC Pack 2008 R2 Service Pack 2</p>
<p class="p18">http://www.microsoft.com/download/en/details.aspx?id=26646</p>
<p class="p13"><br></p>
<p class="p4"><br></p>
<p class="p3"><b>6) install<span class="Apple-converted-space">  </span>inno setup to build installer:<span class="Apple-converted-space"> </span></b></p>
<p class="p12">http://www.xs4all.nl/~mlaan2/ispack/isetup-5.4.0.exe</p>
<p class="p4"><br></p>
<p class="p3"><b>7) GSL for gsl interface<span class="Apple-converted-space">  </span>from</b><span class="Apple-converted-space"> </span></p>
<p class="p5"><span class="s5"><a href="http://sourceforge.net/projects/mingw-cross/files/%5BLIB%5D%20GSL/mingw32-gsl-1.14-1/mingw32-gsl-1.14-1.zip/download">Download Now! mingw32-gsl-1.14-1.zip (3.5 MB)</a></span></p>
<p class="p4"><br></p>
<p class="p3">8) download mercurail for windows from:</p>
<p class="p5"><span class="s5"><a href="http://mercurial.selenic.com/">http://mercurial.selenic.com</a></span></p>
<p class="p4"><br></p>
<p class="p3">9) To download the latest freefem++   tar.gz file contening source form  </p>
    <? pfile($sfile,"Source code") ?>
 <li>  or you can get the latest source from
an anonymous <A href="http://mercurial.selenic.com/" > Mercurial SCM </A> copy with the following unix shell commands&nbsp;:

<p class="p8">hg clone<span class="Apple-converted-space">  </span>http://www.freefem.org/ff++/ff++</p>
<p class="p9">to update do to the last version:</p>
<p class="p8">hg pull</p>
<p class="p8">hg up<span class="Apple-converted-space"> </span></p>
<p class="p4"><br></p>
<p class="p3"><b>To restore, all file build by autoreconf -i command  (automake c):</b></p>
<p class="p7">tar zxvf AutoGeneratedFile.tar.gz </p>

<p class="p3"><b>Finaly, the configure argument are:</b></p>
<p class="p7">cd ff++</p>
<p class="p7">./configure '--enable-download' 'FC=mingw32-gfortran' 'F77=mingw32-gfortran' 'CC=mingw32-gcc' 'CXX=mingw32-g++' '-with-blas=/home/hecht/blas-x86/libgoto2.dll' 'CXXFLAGS=-I/home/hecht/blas-x86' '--enable-generic' '--with-wget=wget' 'MPIRUN=/c/Program Files/Microsoft HPC Pack 2008 R2/Bin/mpiexec.exe'</p>
<p class="p19"><br></p>
<p class="p20">if erreor where building DOC do:</p>
<p class="p7"><span class="s7">t</span>ouch DOC/freefem++doc.pdf</p>
<p class="p4"><br></p>
<p class="p4"><br></p>
<p class="p3">Good Luck …<span class="Apple-converted-space"> </span></p>
<p class="p4"><br></p>
<p class="p4"><br></p>
</body>
</html>
