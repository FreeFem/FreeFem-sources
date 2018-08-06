#!/usr/local/bin/perl
# -----  clean ---
unlink <*.mesh>;
unlink <*.am_fmt>;
unlink <*.mtr>;
unlink <*.bb>;
unlink <*.BB>;
unlink <YG_trace*>;
unlink 'PLOT';

##$f="10 +  1/(1+ 100**(sin(x*3)-y)) ";
#$f1 = "(10*x*x*x+y*y*y) + 10/(1+10**(10*((sin(5*x)-2*y)))) ";
$f1=" sin(3*x)*cos(5*y)+ atan2(0.001,x**2 + y**2 - 0.5 )";
#$f1=" 8*(x-y)**2 + (x+y)**2";

$err=0.1;
$errg=0.01;
$nbiteration=10;
$bamg=$ENV{'bamg'};
$quadoption="";
$bamgoption=" -AbsError -NbJacobi 2  -NbSmooth 5 -anisomax 5 -hmax 0.5   -ratio 2  -nbv 100000 ";
#$bamgoption=" -AbsError -NbJacobi 3  -ratio 2 -anisomax 30";
#$quadoption=" -2q -thetaquad 30 -coef 2";

# ---  change x in $x and y in $y 
$_=$f1;
s/x/\$x/g;
s/y/\$y/g;
$f1="$_;";

print "The function = f1(x,y) = $f1 \n";

#--------------------------
$suffixe=".mesh";
$iteration=0;
$GH="Gh$suffixe";
$TH="Th$iteration$suffixe";
$, = ' ';               # set output field separator
$\ = "\n";              # set output record separator
##  -------------------------------------------------------------
##  --- construction the Geometry file Gh.mesh 

##    8 points on circle of radius r 
open(GH,">$GH")  || die "Can't redirect stdout";
$Pi = 3.14159265358979;

$i8 = 8;
$r = 1;
$c0x = 0;
$c0y = 0;
print GH 'Dimension', 2;
print GH 'MaximalAngleOfCorner 46';

print GH 'Vertices';
print GH   $i8;

# vertex on circle  (counter clock wise)
for ($i = 0; $i < $i8; $i++) {
    $t = $i * $Pi * 2 / $i8;
    print GH $c0x + $r * (cos($t)), $c0y + $r * sin($t), 5;
}

print GH 'Edges',  $i8+1;

print GH 1,5,10;
# edge on circle 
$k =  1;
$j = $i8 - 1;
# previous points
for ($i = 0; $i < $i8; $j = $i++) {
    print GH $k + $j, $k + $i, 5;
}
# previous, current vertex

#   one subdomain, region on left side of the wing
#   because clock wise sens.
print GH 'SubDomain', 2;
print GH 2, 1, 1, 0;
print GH 2, 1, -1, 1;
close GH;
##  -------------- END construct of the geom 
## 

# -- make the DATA file for the mesh to also  save the arguments 
open(BAMG,">DATA_bamg")  || die "Can't open  DATA_bamg";
print BAMG "$quadoption $bamgoption  -g $GH -o $TH  -v 9 -oam_fmt $TH.am_fmt";
close(BAMG);

## constructio the inital  mesh 
!system($bamg) || die "Error in bamg construction of initial  mesh $Th";

##  the adpatation loop 
while ($iteration<$nbiteration) {

    

    $BB="$iteration.BB";   
    $ERRBB="err$iteration.bb";   

    
##  construction of the solution  
    $errsol=0;
    open (TH,"<$TH") ||  die "Can't open  $TH";
    open (BB,">$BB") ||  die "Can't open  $BB";
    open (ERRBB,">$ERRBB") ||  die "Can't open  $ERRBB";


    while (<TH>) {
	if(/^Vertices$/) {
	    $nbv=<TH>;
	    chop($nbv);
	    print BB "2 1 1 $nbv 2";

	    for ($i=1;$i<=$nbv;$i++) {
		($x,$y,$ref)=split(/[\ \t\n]+/, <TH>);
		$f1xy=eval $f1;

		$xx[$i]=$x;
		$yy[$i]=$y;
		$ff[$i]=$f1xy;


		print BB $f1xy ;

	    };
	};
	if(/^Triangles$/) {
	    $nbt=<TH>;
	    chop($nbt);
	    print " Nb of Triangles = $nbt \n";
	    print ERRBB "2 1 $nbt 1";
	    for ($i=1;$i<=$nbt;$i++) {
		($i0,$i1,$i2,$ref)=split(/[\ \t\n]+/, <TH>);

		$x   = ($xx[$i0]+$xx[$i1]+$xx[$i2])/3;
		$y   = ($yy[$i0]+$yy[$i1]+$yy[$i2])/3;
		$fm  = ($ff[$i0]+$ff[$i1]+$ff[$i2])/3;
		$fxy = eval $f1;
		$vv=($fm-$fxy);
		$vv= ($vv<0)?-$vv:$vv;
		print ERRBB $vv;
		$errsol = ($errsol>$vv) ? $errsol : $vv;
	    };
	};

    };
    close TH;
    close BB;
    close ERRBB;
    print " ---------------------------------------------\n";
    print "\n\n Iteration $iteration\n Erreur L_infini = $errsol \n\n";
    print " ---------------------------------------------\n";

##  -----------------------
    
    $MTR="M$iteration.mtr";

    $iteration++;
    $BTH=$TH;
    $TH="Th$iteration$suffixe";

    open(BAMG,">DATA_bamg")  || die "Can't open  DATA_bamg";
    print BAMG "$quadoption  $bamgoption  -MBB $BB -errg $errg -err $err   -b $BTH -o $TH  -v 9 -oM $MTR -oam_fmt $TH.am_fmt -wBB /tmp/tyty ";
    close(BAMG);
    !system($bamg) ||    die "Error in bamg construction of adapted $iteration  mesh $Th";
}

print "Normal End\n";
