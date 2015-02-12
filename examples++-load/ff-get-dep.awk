BEGIN { nl=split(libs,l,"[ ]*"); err= 0;sp=" ";db=0;}
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
	libb=l[i];
	if (match(libb,/^[ ]*$/)) continue; 
	lib=libb;
	nn=sub(/\[/,"",lib);
	nn+=sub(/\]/,"",lib);
	if (ld[lib]=="" && nn ==0) { lerr[err++]=lib ;}
        if(ld[lib]!="" && nn==2) msg=msg "-DWITH_" lib sp;
	if(db) print lib, err;
	msg = msg ld[lib] sp inc[lib] sp; 
    }
    
    if(err) {
	for(i=0; i < err;++i)
	    print " MISSING lib " , lerr[i] 
	print " Check the WHERE-LIBRARY... files ";
	exit 1; 
    }
    else printf("%s\n",msg); 
}