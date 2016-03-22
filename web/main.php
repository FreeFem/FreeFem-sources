<html>
<head>
<link rel="stylesheet" href="ffstyle.css" type="text/css">
<meta http-equiv="Content-type" content="text/html; charset=utf-8">
   <meta name="Description" content="FreeFem++ is a language that allows 
the resolution of partial      differential equation using the finite 
element method">
   <meta name="Keywords" content="language, c++, Multi-platform, free      
software, Navier-Stokes, elasticity, convection-diffusion, heat      
equation, linear elliptic PDE's",MPI,"Scientific computing">

   <title>Freefem++ Home Page (March .  2016)</title>
   <link href="mailto:frederic.hecht@upmc.fr" rev="Author">

</head>

<?php
include 'phpfiles.php'
?>



<body bgcolor="#FFFFFF" text="#000000" link="#333333" vlink="#330066" alink="#330066">
<div class="thetitle"> FreeFem++ v <? echo $fver ?>  <font size=-1> (<? echo $fdate?>) </font> </div> 
<div class="content">
<div class="thema">Introduction</div>
<div class="themaBlog">
</tr></table>

<p>&nbsp;

<table><tr>

<td><img src="images/csSnapSmall.jpg"></td>

<td><b>FreeFem++</b> is a partial differential equation solver.
 It has its own language. freefem scripts can solve multiphysics non linear systems in 2D and 3D.

<p>Problems involving PDE (2d, 3d)  from several branches of physics such as
fluid-structure interactions require interpolations of data on several
meshes and their manipulation within one program. FreeFem++ includes a
fast 2^d-tree-based interpolation algorithm and a language for the
manipulation of data on multiple meshes (as a follow up of 
bamg (now a part of FreeFem++ ).

<p>FreeFem++ is written in C++ and the FreeFem++ language is a C++
idiom.
It runs on Macs, Windows, Unix machines.  FreeFem++ replaces the older <a
href="freefem/fraold.htm">freefem</a> and <a
href="freefem/frap.htm">freefem+</a>.  </td>
</tr></table>

</div> 

If you use <tt size=-2>Freefem++</tt> please cite the following reference in your work (books, articles, reports, etc.):
<a href="http://dx.doi.org/10.1515/jnum-2012-0013" > Hecht, F. New development in FreeFem++. J. Numer. Math. 20 (2012), no. 3-4, 251–265. 65Y15<a> 
</p> </a>
the bibtex is:
<table><tr>
<td>
<blockquote><TT><PRE>
@article {MR3043640,
    AUTHOR = {Hecht, F.}, TITLE = {New development in FreeFem++},
   JOURNAL = {J. Numer. Math.},  FJOURNAL = {Journal of Numerical Mathematics},
    VOLUME = {20}, YEAR = {2012},
    NUMBER = {3-4}, PAGES = {251--265},
      ISSN = {1570-2820}, MRCLASS = {65Y15}, MRNUMBER = {3043640},
}
<PRE><TT/></blockquote></td></tr><table>

<br><br>
<!-- 
       <li><font > <bf>  <font color=black> 
       The 6th   tutorial and Workshop on FreeFem++ , 
 held december 9th,10th and 11th, 2014, in Paris at Universit&eacute Pierre
et Marie Curie, Barre 16-15, 3ieme,  4 place Jussieu, Paris</bf>  
<A HREF="https://www.ljll.math.upmc.fr/~hecht/ftp/ff++days/2014/Schedule.html"> Schedule, All Presentation, examples </A> and
<A HREF="https://www.ljll.math.upmc.fr/~hecht/ftp/ff++days/2014/.">  directory this all data  </A>.
 </font></li> 

           <li><font > <bf>  <font color=black> 
       The 7th   tutorial and Workshop on FreeFem++ 
 <A HREF="https://www.ljll.math.upmc.fr/FreeFem++/">  (inscription here)</A>, 
 held december 15th and 16th, 2015, in Paris at Universit&eacute Pierre
et Marie Curie, Barre 16-15, 3ieme,  4 place Jussieu, Paris</bf>  
<A HREF="https://www.ljll.math.upmc.fr/~hecht/ftp/ff++days/2015/Schedule.html"> Schedule, All Presentation, examples </A> and
<A HREF="https://www.ljll.math.upmc.fr/~hecht/ftp/ff++days/2015/.">  directory this all data  </A>.
 </font></li> 
  
    
  </ul>
  
</div> --> 
<div class="thema"> HPC and FreeFem++  </div>
<ul>
<li> Solver a problem with  one  billion of unknows  in  2 minutes! The computation will done in 2 minutes on the  Curie Thin Node@CEA machine (6144 coeurs de 4Go de mémoire chacun),  Thank to P. Jolivet, and F. Nataf. </li>
<li>  In construction </li>


</ul>

<div class="thema">Some FreeFem++ presentation (with useful information):
</div>
<ul>
    <li> <font > <bf>  <font color=blue>  Graduate Course: </bf> </font>
<A HREF="https://www.fields.utoronto.ca/programs/scientific/15-16/scientificcomputing/graduate/">  
 An introduction to scientific computing using free software FreeFem++ <A/> , 
 <A HREF="https://www.fields.utoronto.ca/"> The Fields Institute for 
Research in Mathematical Sciences <A/>,  Toronto, Canada , 7-17 March . 2016   
<A HREF="https://www.ljll.math.upmc.fr/~hecht/ftp/ff++/2016-Fields"> The directory with all data  </A>,  
<A HREF="https://www.ljll.math.upmc.fr/~hecht/ftp/ff++/2016-Fields.zip"> The zip of the directory with all data  </A> 
 </li>

 
         <li><font > <bf>  <font color=black> <font color=black> 
       The 7th   tutorial and Workshop on FreeFem++ </font>
 <A HREF="http://www.ljll.math.upmc.fr/FreeFem++">  (inscription here)</A>, 
 held december 15th and 16th, 2015, in Paris at Universit&eacute Pierre
et Marie Curie, Barre 16-15, 3ieme,  4 place Jussieu, Paris</bf>  
<A HREF="http://www.ljll.math.upmc.fr/~hecht/ftp/ff++days/2015/Schedule.html"> Schedule, All Presentation, examples </A> and
<A HREF="http://www.ljll.math.upmc.fr/~hecht/ftp/ff++days/2015/.">  directory this all data  </A>.
 </font></li> 

   <li>  FreeFem++, 
<A HREF="http://www.math.iitb.ac.in/~neela/CIMPA/home.html">   Cimpa Summer School on Current Research in FEM  at IIT Bombay, India <A/> 6-17 July . 2015   
<A HREF="https://www.ljll.math.upmc.fr/~hecht/ftp/ff++/2015-cimpa-IIT"> The directory with all data  </A> 
 </li>

  <li> 

 Mini Cours FreeFem++ <A HREF="http://www.maths.ox.ac.uk"> Maths departement, Universty of  Oxford, England  <A/> 16-20 March . 2015   
<A HREF="http://www.ljll.math.upmc.fr/~hecht/ftp/ff++/2015-oxford"> The directory with all data  </A> 
 </li>
       <li><font > <bf>  <font color=black> 
       The Sixth   tutorial and Workshop on FreeFem++ 
 <A HREF="http://www.ljll.math.upmc.fr/FreeFem++"> 
 (inscription here)</A>, 
 held december 9th,10th and 11th, 2014, in Paris at Universit&eacute Pierre
et Marie Curie, Barre 16-15, 3ieme,  4 place Jussieu, Paris</bf>  
<A HREF="http://www.ljll.math.upmc.fr/~hecht/ftp/ff++days/2014/Schedule.html"> Schedule, All Presentation, examples </A> and
<A HREF="http://www.ljll.math.upmc.fr/~hecht/ftp/ff++days/2014/.">  directory this all data  </A>.
 </font></li> 

<li>  FreeFem++, 
<A HREF="http://cadiz-numerica-2013.uca.es">  Cadiz numerica 2013 <A/>
University of Cadiz, Spain ,  26-28 Juin . 2013   <A HREF="http://www.ljll.math.upmc.fr/~hecht/ftp/FF-conf/cadiz-numerica-2013"> All presentation and script  </A> </li>

<li>  FreeFem++, University of Houston, TX, USA,  8,-9, feb. 2013   <A HREF="http://www.ljll.math.upmc.fr/~hecht/ftp/FF-conf/UoH-2013"> All presentation and script  </A> </li>
   
 <li>     Cours de presentation FreeFem++ 
    <A HREF=http://master-modsim.univ-rennes1.fr/"> Master I,  Modélisation et calcul scientifique, Rennes I  </A>:
    <A HREF="ftp/Cours/Rennes"> Site du cours  </A>.</li>     
         <li>
 A seminar at <A HREF="http:http://departamento.us.es/edan/">  Departamento de Ecuaciones Diferenciales y Análisis Numérico
 </A> de la Universidad de Sevilla, Espagne, Jan, 2011   <A HREF="http://www.ljll.math.upmc.fr/~hecht/ftp/FF-conf/FF-Seville-2011-MPI-Schwarz-DDM.zip"> Coarse grid Schwarz DDM parellel solver in FreeFem++ (y &epsilon;) </A> </li>

  
        <li>
 A seminar at <A HREF="http://www-math.univ-paris13.fr/laga/">  laga  </A>
 at Villetaneuse, france , jan, 2010   <A HREF="http://www.ljll.math.upmc.fr/~hecht/ftp/FF-conf/ff++-laga-2010.pdf">FreeFem++, a tool to solve PDE’s numerically. </A>
</li>
 
      <li>
 A presentation at <A HREF="http://lma.univ-pau.fr/meet/mamern09/">  MAMERN 09 .</A>
 at Pau, france , june 11th, 2009,   <A HREF="http://www.ljll.math.upmc.fr/~hecht/ftp/FH-Mamerm09.pdf">the slide  Error indicator and mesh
adaption, in FreeFem++ </A>.
</li>
      <li>
 A presentation and   all the data files at <A HREF="http://math.tkk.fi/numericsyear/fefair/"> Finite element fair .</A>
Helsinki University of Technology 
Institute of Mathematics, june 5-6Th 2009 , the zip file :<A HREF="http://www.ljll.math.upmc.fr/~hecht/ftp/HECHT-FEFE09-dir.zip "> FreeFem++, 3d tools for 
PDE simulation </A>.</li>
</li>
   
     <li>
 <A HREF="http://www.cimpa-icpam.org/index.php"> CIMPA</A>-UNESCO-
 <A HREF="http://www.univ-ag.fr/aoc/?guadeloupe09">GUADELOUPE School  </A>
January, 03-18, 2009, Pointe-à-Pitre, the <A HREF="ftp/CIMPA/CIMPA-Guadeloupe-FF.pdf"> pdf  of my lecture (742Ko)</A> on FreeFem++ and the  <A HREF="ftp/CIMPA/FF-CIMPA.zip">  archive (.zip) of all the examples (6.2Mo) </A>
</li>
 
       <li>  Numerical modeling of Geophysical Flows by Finite Element techniques with FreeFem++ <A HREF="http://institucional.us.es/doc-course-imus/Unit2.html" >  IMUS 2010 Univerty of Seville</A>,    Spain, <A HREF="ftp/IMUS-Seville"> Site of the cours  </A>  
      . </li>

         <li>
 My seminar at <A HREF="http:http://departamento.us.es/edan/">  Departamento de Ecuaciones Diferenciales y Análisis Numérico
 </A> de la Universidad de Sevilla, Espagne,  , janvier, 2011   <A HREF="http://www.ljll.math.upmc.fr/~hecht/ftp/FF-conf/FF-Seville-2011-MPI-Schwarz-DDM.zip"> Coarse grid Schwarz DDM parellel solver in FreeFem++ (y &epsilon;) </A> (to install
  my  FreeFem++ mpi version on your Mac 10.6 (Snow leopard)   <A HREF="http://www.ljll.math.upmc.fr/~hecht/ftp/FF-conf/InstallMac10.6.html" > See this page </A>). 
</li>

</ul> 

<div class="thema">Related software:</div>
</div>
<ul>
<li>     Emc2 (Editeur de Maillage 2d) <A HREF="http://www.ljll.math.upmc.fr/~hecht/ftp/emc2/"> la dernière version (LJLL) </A> et  <A HREF="http://www.ljll.math.upmc.fr/~hecht/ftp/emc2/RTemc2_gb.pdf"> english documentation in  pdf </A>. 
<!--    et  sur le site de l'INRIA:  <A HREF="http://www-c.inria.fr/gamma/cdrom/www/emc2/fra.htm">  à l'inria en français </A>,
    <A HREF="http://www-c.inria.fr/gamma/cdrom/www/emc2/eng.htm">  in english </A>. -->
     </li> 
</ul>


<br><br>
<div class="thema">Examples 2d </div>
<div class="themaBlog">  
 <ul>

<li><a href="freefem/ff.mp4">A small movie (340Kb)</a>&nbsp;:
Cool air (green) comes from the lower left and mix with hot air
(magenta), the right boundary is free. This is
Navier-Stokes-Boussinesq integrated with P1-bubble P1 mixte finite
element.

</ul>

<ul>
<li>A very small example 2d of how to solve the Poisson equation on a L
shape&nbsp;:

<table><tr>

<td>
<blockquote><TT><PRE>
border aaa(t=0,1){x=t;y=0;};
border bbb(t=0,0.5){x=1;y=t;};
border ccc(t=0,0.5){x=1-t;y=0.5;};
border ddd(t=0.5,1){x=0.5;y=t;};
border eee(t=0.5,1){x=1-t;y=1;};
border fff(t=0,1){x=0;y=1-t;};
mesh Th = buildmesh (aaa(6) + bbb(4) + ccc(4) +ddd(4) + eee(4) + fff(6));
fespace Vh(Th,P1);  <font color=red> //  to change P1 in P2 to make P2 finite element.</font>
Vh u=0,v;
func f= 1;
func g= 0;
int i=0;
real error=0.1, coef= 0.1^(1./5.);
problem Probem1(u,v,solver=CG,eps=-1.0e-6) =
    int2d(Th)(  dx(u)*dx(v) + dy(u)*dy(v)) 
  + int2d(Th) ( v*f ) 
  + on(aaa,bbb,ccc,ddd,eee,fff,u=g)  ;
  
for (i=0;i< 10;i++)
{   
  real d = clock();
  Probem1; <font color=red>//  solve the problem </font>
  plot(u,Th,wait=1);
  Th=adaptmesh(Th,u,inquire=1,err=error);
  error = error * coef;
} ;

</PRE></TT></blockquote>
</td>

<td>
<img SRC="images/u4.jpg" BORDER=0 height=250 width=250><img
SRC="images/th4.jpg" BORDER=0 height=250 width=250> </p>
<p align=CENTER> Solution on adapted mesh and associated mesh. </p>
</td>

</tr>
</table>

</ul>

</div> 


<div class="thema">Examples 3d </div>
<div class="themaBlog">  

<ul>
<li>A very small example of how to solve the Stokes equation 3d on cube
shape&nbsp;:

<table><tr>

<td>
<blockquote><TT><PRE><font color=black>
load "msh3" load "medit" <font color=red> // dynamics load tools for 3d.</font>
int nn=8;
mesh Th2=square(nn,nn);
fespace Vh2(Th2,P2);  Vh2 ux,uz,p2;
int[int] rup=[0,2],  rdown=[0,1], rmid=[1,1,2,1,3,1,4,1];
real zmin=0,zmax=1;
mesh3 Th=buildlayers(Th2,nn,
  zbound=[zmin,zmax],  reffacemid=rmid, 
  reffaceup = rup,     reffacelow = rdown);
  
medit("c10x10x10",Th); <font color=red> // see the 3d mesh with medit software</font>
fespace VVh(Th,[P2,P2,P2,P1]);

<font color=blue> macro Grad(u) [dx(u),dy(u),dz(u)] // EOM</font>
<font color=blue> macro div(u1,u2,u3) (dx(u1)+dy(u2)+dz(u3))  //EOM</font>

VVh [u1,u2,u3,p];
VVh [v1,v2,v3,q];
  
solve vStokes([u1,u2,u3,p],[v1,v2,v3,q]) = 
  int3d(Th,qforder=3)( Grad(u1)'*Grad(v1) +  Grad(u2)'*Grad(v2) +  Grad(u3)'*Grad(v3)
                  - div(u1,u2,u3)*q - div(v1,v2,v3)*p + 1e-10*q*p ) 
  + on(2,u1=1.,u2=0,u3=0) + on(1,u1=0,u2=0,u3=0) ;
 plot(p,wait=1, nbiso=5); <font color=red> // a 3d plot of iso  pressure. in progress... march 2009</font>
<font color=red> //  to see the 10 cut plan in 2d </font>
for(int i=1;i<10;i++)
{
 real yy=i/10.; <font color=red>// compute yy.</font>
 <font color=red> // do 3d -> 2d interpolation.</font>
 ux= u1(x,yy,y); uz= u3(x,yy,y);  p2= p(x,yy,y);
 plot([ux,uz],p2,cmm=" cut y = "+yy,wait= 1);
}
</PRE></TT></blockquote>
</td>

<td>
<img SRC="images/Stokes3d.jpg" BORDER=0 height=250 width=250><img
SRC="images/Th-Stokes3d.jpg" BORDER=0 height=250 width=250> </p>
<p align=CENTER> Solution on cup plan y=0.5 and mesh 10x10x10  and associated mesh. </p>
</td>

</tr>
</table>

</ul>
<!-- ---------------------------------------------------------------------- -->
<div class="thema"> Ongoing Work </div>
<div class="themaBlog">  

<ul>
<li>
FreeVol: Finite Vol technics in FreeFem++ for hyperbolic PDEs </li>
<li>
3D implementation: new solver, new mesh tools, new kind of finite
element </li>

<li> Stabilize and test all parallel linear solver interface</li>

</ul>

</div>

<!-- ---------------------------------------------------------------------- -->

<!-- ---------------------------------------------------------------------- -->
</div>   
<div class="thema"> Download,  The current version of <i>FreeFem++</i> is <? echo $fver ?> </div>
<div class="themaBlog">  


<ul> <li>  You can get the latest source from
an anonymous <A href="http://mercurial.selenic.com/" > Mercurial SCM </A> copy with the following unix shell commands&nbsp;:
<p><TT><font color=black size=+0>
hg  clone  http://www.freefem.org/ff++/ff++
</font></tt></p>
</li> 
Self-contained archives for all other systems&nbsp;:
<P>

<? download($adown) ?>
</p>
</menu>

<li>   <p>This <a href="ftp"> directory</A> contains all the different
versions of <i>FreeFem++</i>.</font> </li> 
</ul>
</div>    
 
<address></address>
<!-- hhmts start --> Last modified: 2 juin  2014 <!-- hhmts end -->
</body> </html>
</body>
</html>
