#!/usr/bin/perl
# -- so option
 
##$f="10 +  1/(1+ 100**(sin(x*3)-y)) ";
$f = "10 + sin(x/10)*cos(y/3)";
$err=0.001;
$errg=0.05;
$nbiteration=3;
$bamg="../../src/bamg/bamg";
$quadoption="";
$bamgoption=" -AbsError -NbJacobi 3  -ratio 2 -anisomax 30";
$quadoption=" -2q -thetaquad 30 -coef 2";

# ---  change x in $x and y in $y 
$_=$f;
s/x/\$x/g;
s/y/\$y/g;
$f="$_;";
print "The function = f(x,y) = $f \n";

#--------------------------
$suffixe=".mesh";
$iteration=0;
$GH="Gh$suffixe";
$TH="Th$iteration$suffixe";
$, = ' ';               # set output field separator
$\ = "\n";              # set output record separator
##  -------------------------------------------------------------
##  --- construction the Geometry file Gh.mesh 
##  a naca0012 wing in a big circle radius 5 
##   20 points on the up wing 
##   20 points on the down wing 
##    8 points on circle of radius r 
open(GH,">$GH")  || die "Can't redirect stdout";
$Pi = 3.14159265358979;
$i20 = 20;
$i8 = 8;
$r = 5;
$c0x = 1;
$c0y = 0;
print GH 'Dimension', 2;
print GH 'MaximalAngleOfCorner 46';

print GH 'Vertices';
print GH $i20 + $i20 + $i8;

# the vertex on naca012  wing ( clock wise)
for ($i = -$i20; $i < $i20; $i++) {
    $X = $i / $i20;
    $X = $X ** 4;

    $t = 1.008930411365 * $X;
    $Y = 5 * .12 * (0.2969 * sqrt($t) - 0.126 * $t - 0.3516 * $t ** 2 + 0.2843

      * $t ** 3 - 0.1015 * $t ** 4);
    if ($i < 0) {
        $Y = -$Y;
    }
    print GH $X, $Y, 3;
}

# vertex on circle  (counter clock wise)
for ($i = 0; $i < $i8; $i++) {
    $t = $i * $Pi * 2 / $i8;
    print GH $c0x + $r * (cos($t)), $c0y + $r * sin($t), 5;
}

print GH 'Edges', $i20 + $i20 + $i8;

#  edge on wing  
$k = 1;
$j = $i20 + $i20 - 1;
# previous points
for ($i = 0; $i < $i20 + $i20; $j = $i++) {
    print GH $j + $k, $i + $k, 3;
}
# previous, current vertex

# edge on circle 
$k = $i20 + $i20 + 1;
$j = $i8 - 1;
# previous points
for ($i = 0; $i < $i8; $j = $i++) {
    print GH $k + $j, $k + $i, 5;
}
# previous, current vertex

#   one subdomain, region on left side of the wing
#   because clock wise sens.
print GH 'SubDomain', 1;
print GH 2, 1, 1, 0;
close GH;
##  -------------- END construct of the geom 
## 

# -- make the DATA file for the mesh to also  save the arguments 
open(BAMG,">DATA_bamg")  || die "Can't open  DATA_bamg";
print BAMG "$quadoption $bamgoption  -g $GH -o $TH  -v 9";
close(BAMG);

## constructio the inital  mesh 
!system($bamg) || die "Error in bamg construction of initial  mesh $Th";

##  the adpatation loop 
while ($iteration<$nbiteration) {

    
    $BB="$iteration.bb";
    
##  construction of the solution  
    $errsol=0;
    open (TH,"<$TH") ||  die "Can't open  $TH";
    open (BB,">$BB") ||  die "Can't open  $BB";
    open (PLOT,">PLOT") || die "Can't open PLOT";
    while (<TH>) {
	if(/^Vertices$/) {
	    $nbv=<TH>;
	    chop($nbv);
	    print BB "2 1 $nbv 2";
	    for ($i=1;$i<=$nbv;$i++) {
		($x,$y,$ref)=split(/[\ \t\n]+/, <TH>);
		$fxy=eval $f;

		$xx[$i]=$x;
		$yy[$i]=$y;
		$ff[$i]=$fxy;

		print BB $fxy;
	    };
	};
	if(/^Triangles$/) {
	    print " Nb of Triangles = $nbt \n";
	    $nbt=<TH>;
	    chop($nbt);
	    for ($i=1;$i<=$nbt;$i++) {
		($i0,$i1,$i2,$ref)=split(/[\ \t\n]+/, <TH>);
		print PLOT "$xx[$i0] $yy[$i0] $ff[$i0]";
		print PLOT "$xx[$i1] $yy[$i1] $ff[$i1]";
		print PLOT "$xx[$i2] $yy[$i2] $ff[$i2]";
		print PLOT "$xx[$i0] $yy[$i0] $ff[$i0]";
		print PLOT "";

		$x   = ($xx[$i0]+$xx[$i1]+$xx[$i2])/3;
		$y   = ($yy[$i0]+$yy[$i1]+$yy[$i2])/3;
		$fm  = ($ff[$i0]+$ff[$i1]+$ff[$i2])/3;
		$fxy = eval $f;
		$vv=($fm-$fxy);
		$vv= ($vv<0)?-$vv:$vv;
	#	print " $i0 $i1 $i2 $xx[$i0] $xx[$i1] $xx[$i2] ";
	#	print " $i $x $y $fm $fxy err= $errsol diff=$vv";
		$errsol = ($errsol>$vv) ? $errsol : $vv;
	    };
	};

    };
    close TH;
    close BB;
    close PLOT;
    print " ---------------------------------------------\n";
    print "\n\n Iteration $iteration\n Erreur L_infini = $errsol \n\n";
    print " ---------------------------------------------\n";

##  -----------------------
    
    $iteration++;
    $BTH=$TH;
    $TH="Th$iteration$suffixe";
    
    open(BAMG,">DATA_bamg")  || die "Can't open  DATA_bamg";
    print BAMG "$quadoption  $bamgoption  -Mbb $BB -errg $errg -err $err   -b $BTH -o $TH  -v 9";
    close(BAMG);
    !system($bamg) ||    die "Error in bamg construction of adapted $iteration  mesh $Th";
}

print "Normal End\n";
