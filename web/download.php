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

   <title>Freefem++ download  Page (May.  2014)</title>
   <link href="mailto:frederic.hecht@upmc.fr" rev="Author">

</head>

<?php
include 'phpfiles.php'
?>



<body bgcolor="#FFFFFF" text="#000000" link="#333333" vlink="#330066" alink="#330066">
<div class="thetitle"> FreeFem++ v <? echo $fver ?>  <font size=-1> (<? echo $fdate?>) </font> </div> 
<div class="content">
<div class="thema"> Download,  The current version of <i>FreeFem++</i> is <? echo $fver ?> </div>
<div class="themaBlog">  


<ul> <li>  You can get the latest source from
an anonymous <A href="http://mercurial.selenic.com/" > Mercurial SCM </A> copy with the following unix shell commands&nbsp;:
<p><TT><font color=black size=+0>
hg  clone  http://www.freefem.org/ff++/ff++
</font></tt></p>
</li> 
<li>     <? pfile($sfile,"Source code") ?> </li> 
</div>

<div class="thema"> Download precompile version  </div>
<div class="themaBlog">  
    <? pfile($w64file,"Windows 64bit (in test)") ?>
    <? pfile($wfile,"Windows 32bit") ?>
    
	<? pfile($m11file,"MacOS 10.11") ?>
	<? pfile($m10file,"MacOS 10.10") ?>
	<? pfile($m9file,"MacOS 10.9") ?>
	<? pfile($m8file,"MacOS 10.8") ?>
	<? pfile($m7file,"MacOS 10.7") ?>	
	<? pfile($m6file,"MacOS 10.6") ?>	
	<? pfile($mfile,"MacOS 10.4") ?>	
	
<p></p>
<div class="thema"> Self-contained archives for all other systems&nbsp;: </div>
<div class="themaBlog">  
<? download($adown) ?>
</p>
</menu>

<div class="thema"> All  the  versions of <i>FreeFem++</i></div>
<div class="themaBlog">  

<li>   <p> are <a href="ftp"> in this place</A>. </li> 
</ul>
</div>    
 
<address></address>
<!-- hhmts start --> Last modified: 30 May 2014 <!-- hhmts end -->
</body> </html>
</body>
</html>
