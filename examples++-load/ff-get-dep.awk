BEGIN { nl=split(libs,l," *"); err= 0;sp=" ";db=0;}
$2=="LD" { 
    if( ld[$1]=="" )  {
	a=$1;
	sub(/.* LD /,"");
	ld[a]=$0 sp; 
	if(db) print a,$0
    }}
$2=="INCLUDE" { if( inc[$1]=="" ) {
	a=$1;sub(/.* INCLUDE /,"");
	inc[a]=$0 sp;
	if(db)	print a,$0;
}}
END {
    msg=cpp sp;
    for(i=1; i<=nl;++i)
    {
	libbb=l[i];
        m=split(libbb,ll,"[]|[]");
        first = 0;
	for(j=1; j<=m; ++j)
	{
	    msg1="";
	    lib = ll[j];
	    nn=1;
	    if(libbb==lib) nn=0;
	    if (ld[lib]=="" && nn ==0) { lerr[err++]=lib ;}            
	    if(ld[lib]!="" && nn!=0)
	    {
		if( first==0) {
		    msg1="-DWITH_" lib sp ld[lib] sp inc[lib] sp;
		}
		first ++;
	    }
	    if(ld[lib]!="" && nn==0)
		msg1 =  ld[lib] sp inc[lib] sp;
	    msg =msg msg1; 
	    if(db) print " ###" libbb,  lib, err, first ,j, m, msg1;
	}

 
    }
    
    if(err) {
	printf  "\t\t  MISSING lib "> "/dev/stderr" 
	for(i=0; i < err;++i)
	    printf "%s, ", lerr[i]> "/dev/stderr"; 
	print "         Check the WHERE-LIBRARYfiles " > "/dev/stderr";
	exit 1; 
    }
    else printf("%s\n",msg); 
}