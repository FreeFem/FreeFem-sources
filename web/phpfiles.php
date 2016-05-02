<?
$ver=0;$ffile="";$fver="";$fdate="";
$wver=0;$wfile="";$wwver="";
$w64ver=0;$w64file="";$w64wver="";
$mver=0;$mfile="";$mmver="";
$m6ver=0;$m6file="";$m6mver="";
$m7ver=0;$m7file="";$m7mver="";
$m8ver=0;$m8file="";$m8mver="";
$m9ver=0;$m9file="";$m9mver="";
$m10ver=0;$m10file="";$m10mver="";
$m11ver=0;$m11file="";$m11mver="";
$sver=0;$sfile="";$ssver="";
$lstfiles="";
$adown = array();


function sizefile($f) 
{
	$file="ftp/".$f;
	$z=round(filesize($file)/1024/1024,1);
	echo " ( ".  $z ." Mb,  ".date ("M d, Y H:i:s.", filemtime($file)). ") \n" ;
    
}
function ahref($f,$txt) 
{
	$file="ftp/".$f;
	echo "<a href=\"".$file."\">$txt</a>";   
}
function print_down($item, $key)
{
    ahref($item,$item) ;
    echo " &nbsp;&nbsp;".sizefile($item) . " <br />\n";
}
function ppfile($f,$ff) 
{
   if (strlen($f) >0) 
   { 
  echo  "<bf>".$ff.":</bf> " ;
  ahref($f,$f) ;
    echo " &nbsp;&nbsp;".sizefile($f) ;

   }
  else
  echo  $ff . ":  do not exist ??? ";
}

function pfile($f,$ff) 
{
   if (strlen($f) >0) 
   { 
  echo  "<li> <bf>".$ff.":</bf> " ;
  ahref($f,$f) ;
    echo " &nbsp;&nbsp;".sizefile($f) . " <br />\n";
    echo  "</li>";
   }
  else
  echo "<li>". $ff . ":  do not exist ???</li> ";
}
function download($adown) 
{
	array_walk($adown,'print_down');
}


if ($handle = opendir('ftp/')) {
    while (false !== ($file = readdir($handle))) {
        if ($file != "." && $file != ".." && !is_dir($file)) {
        	if ( preg_match( "/[v-]([1-9])\.([0-9]+)((-[0-9]+)?)/", $file, $matches ) )
        	{
        		$lver = $matches[1] + $matches[2]*0.001 - $matches[3]*0.000001 ;
        	  //          echo "$file $matches[0] -- $matches[1] -- $matches[3] -- $matches[3].      $ver    ...   ";
          //  echo  (filesize($file)/1024/1024)." Mb,  ".date ("F d Y H:i:s.", filemtime($file)). " " ;
           // echo "<p>";
            if($ver < $lver)
            {
            	$ver=$lver;
            	$ffile=$file;
            	$fdate=date ("F d Y H:i:s.", filemtime("ftp/".$file));
            	$fver = $matches[1].".".$matches[2].$matches[3]; 
            } 
            if(preg_match( "/\.exe/", $file, $mm )  && ($wver < $lver || $w64ver < $lver)  ) 
            {
                if(preg_match( "/win64/", $file, $mm ) )
                {
            	$w64file=$file;
            	$w64ver=$lver;
                }
                else
                {
            	$wfile=$file;
            	$wver=$lver;
                    
                }
             	
        	}
            if(preg_match( "/_10\..*6.*\.pkg/", $file, $mm )  && $m6ver < $lver) 
            {
            	$m6file=$file;
            	$m6ver=$lver;
             	
        	}
            if(preg_match( "/_10\..*7.*\.pkg/", $file, $mm )  && $m7ver < $lver) 
            {
            	$m7file=$file;
            	$m7ver=$lver;
             	
        	}
            if(preg_match( "/_10\..*8.*\.pkg/", $file, $mm )  && $m8ver < $lver) 
            {
            	$m8file=$file;
            	$m8ver=$lver;
             	
        	}
            if(preg_match( "/_10\..*9.*\.pkg/", $file, $mm )  && $m9ver < $lver) 
            {
            	$m9file=$file;
            	$m9ver=$lver;
             	
        	}
             if(preg_match( "/_10\..*10.*\.pkg/", $file, $mm )  && $m10ver < $lver) 
            {
            	$m10file=$file;
            	$m10ver=$lver;
             	
        	}
            if(preg_match( "/_10\..*11.*\.pkg/", $file, $mm )  && $m11ver < $lver) 
            {
            	$m11file=$file;
            	$m11ver=$lver;
             	
        	}
        	
           if(preg_match( "/niversal*\.pkg/", $file, $mm )  && $mver < $lver) 
            {
            	$mfile=$file;
            	$mver=$lver;
             	
        	}
       	
            if(preg_match("/freefem\+\+-[1-9]\.[0-9]+(-[0-9]+)?\.tar\.gz/", $file, $mm )  && $sver < $lver) 
            {
            	$sfile=$file;
            	$sver=$lver;
             	
        	}
 
           
            if(preg_match( "/MacOsX/", $file, $mm )  && $wver < $lver) 
            {
            	$mofile=$file;
             	
        	}
        	preg_match( "/^.*[v-][1-9]\.[0-9]+[-0-9]*(.*)$/", $file, $matche2 );
        	$kf = $matche2[1];
        	if (array_key_exists($kf, $adown))
        	  {
        	  	preg_match( "/[v-]([1-9])\.([0-9]+)([-0-9]*)/",$adown[$kf] , $matche3 );
        	  	$mlver = $matche3[1] + $matche3[2]*0.001 - $matche3[3]*0.000001 ; 
        	  	if( $mlver < $lver )
        	  	 {  //echo " push $kf -> $file ( $lver ) </br>\n ";
        	  	 	$adown[$kf] = $file;}
        	  }
        	else 
        	{  //echo " new push $kf -> $file  ( $lver ) </br>\n";
        	  $adown[$kf] = $file;
        	}
        }}
    }
    closedir($handle);
    rsort($adown);
   // //echo " </br> </br> ************* </br> </br> ";
    //print_r($adown); 
    }
  
    // echo " .... $ver $ffile $fver <p>" ;
?>

