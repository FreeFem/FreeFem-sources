#!/usr/bin/perl 
sub FileInfo {
    my $file = shift;
   my $size = ( -s $file );
   my $size= int(($file/1024)*10)/10;
   my $date_str= localtime((stat($file))[9]);    
   return  "$size Kb $date_str";
}
use File::Basename;
$imagedir="./icons";
$dir=" il manque un parametre $#ARGV ";
if ( $#ARGV == 0 )  { 
 $dir=$ARGV[0];
}
( -d $dir ) ||  die " Erreur $dir n est pas un directory ! "; 
open(file,"ls -1 $dir/[Ff]ree*++*|") ||  die " Erreur ouverture  ";
$kk=0;
$k=0;
$VV = "00000 00000 00000";
FILE:
while (<file>) {
    chop;
   # next FILE if  /^l/  ; # un  lien  on saut
   # next FILE if  /^d/  ; # un  dir  on saut
    next if ( -d $_ );
    next if ( -l $_ );
 
#   ($mode,$nhdcopi,$name,$projet,$size,$mois,$jour,$an,$fichier)=split;
    $fichier=$_; 
    $dir=dirname($fichier);
    $file=basename($fichier);
    $f=$file;
    $p=$file;
    $v=$file;
    $s=$file;
#    $v =~ s/.*\+\+[.v-]?//;
    $v =~ s/^(.*\+\+[.v-]?)([1-9][-.0-9]*[0-9])(.*)$/\2/;
    $p =~ s/^(.*\+\+[.v-]?)([1-9][-.0-9]*[0-9])(.*)$/\1/;
    $s =~ s/^(.*\+\+[.v-]?)([1-9][-.0-9]*[0-9])(.*)$/\3/;
    ($va,$vb,$vc) = split(/[-.]/,$v);
    $clef="$p $s";
    # print " -- $p $v $s $clef  \n";
    if (  $iclef{$clef} == "" ) { $iclef{$clef}=$kk++;}
    $i=$iclef{$clef};
    $pi[$i]=$p;
    $si[$i]=$s;
    $vi[$i] = sprintf("%05d %05d %05d;%s",$va,$vb,$vc,$vi[$i]);
    $VI=  sprintf("%05d %05d %05d",$va,$vb,$vc);
    if($VV le $VI) { $VV=$VI;} 
    $fk[$k] = $_; 
    $ki[$i] = "$k;$ki[$i]";
    $k++;

#    @item=split(/-/);
}
($vv1,$vv2,$vv3) = split(/ /,$VV);
$ver = sprintf("%d.%d-%d",$vv1,$vv2,$vv3);
print " last version: $ver \n";
for ($i=0;$i<$kk;$i++)
{
  # print " $i --- \n" ;
  # print "  $pi[$i]\n";
  # print "  $si[$i]\n";
  # print "  $vi[$i]\n";
   @vvv=split(/;/,$vi[$i]);
   @kkk=split(/;/,$ki[$i]);
    $jj=0;
    for ($j=0;$j<$#vvv;$j++)
    {
	if ( $vvv[$jj] le $vvv[$j])
	{ $jj=$j;}
    }

 #   print " $vvv[$jj]  $kkk[$jj] ??? \n";
    $fff=$fk[$kkk[$jj]];
#   print "$fff \n";
   $filei[$i] = $fff; 
   $tt[$i] = $si[$i] ;
}

# add munal file
for ("HISTORY","manual-full.pdf","manual.pdf","")
{ 
 # print "$dir/$_ ", -f "$dir/$_", "\n";
  next if !( -f "$dir/$_") ;
  $filei[$kk]="$dir/$_";
  $tt[$kk]="$_";
  $kk++; 
}

#  clean the tt array
for (@tt)
{
	s/^.exe$/Setup file for Window 95,98,NT,2000, XP/;
	s/^_MacOsX.tgz$/MacOs 10.3  Cocoa OpenGL (tar+gzip)/;
	s/^_MacOS.sit$/MacOS 9 Powerpc \(carbon StuffIt Archive\)/;
	s/^.tar.gz$/All the sources  (tar+gzip)/;
	s/[.]tar[.]gz$/ (tar+ gzip file)/;
	s/[.]tgz$/ (tar+ gzip file)/;
	s/[.]zip$/ (zip file)/;
	s/[.]pdf$/ (pdf file)/;
	s/^manual-full/Full new manual /;
	s/^manual /Old manual /;
        s/_/  /g;
}
# ici generation de la liste des telechargements
#  $tt[i]    : file caracteristique 
#  $filei[i] : file path 
$download = ""; 
for ($i=0;$i<=$#filei;$i++)
{
    $file = $filei[$i];
#    ($mode,$nhdcopi,$name,$projet,$size,$mois,$jour,$an,$fichier)=split;
    $dir=dirname($file);
    $bfile=basename($file);
    $size = ( -s $file );
    $size= int(($size/1024)*10)/10;
    $date_str= localtime((stat($file))[9]);    
#    $now_string = strftime "%a %b %e %H:%M:%S %Y", localtime ;

   # $date_str = strftime "%a %b %e %H:%M:%S %Y", 
   #$write_secs[0],  $write_secs[1],  $write_secs[2],  $write_secs[4],  $write_secs[5]   ;
#    print "\n @write_secs $date_str \n";
#    $date_str = ctime((stat($file))[9]);
#    $date_str =scalar localtime($write_secs);    
   $download .=  "\n<li>  $tt[$i] \n   <a href=\"$file\"> $bfile </a>\n  <font size=-1>  $size Kb $date_str  </font> </li> \n";
}
$manpdf=FileInfo($dir."/manual-full.pdf");
$omanpdf=FileInfo($dir."/manual.pdf");
$dateofday=localtime();
#print " $download " ; 
open(fout,">freefem++.htm");
while (<STDIN>)
{
    s/\@ver\@/$version/;
    s/\@imagedir\@/$imagedir/;
    s/\@dateofday\@/$dateofday/;
    s/\@download\@/$download/;
    s/\@dir\@/$dir/;
    s/\@manpdf\@/$manpdf/;
    s/\@omanpdf\@/$omanpdf/;
    print fout $_; 
}
close(file);
close(fout);

