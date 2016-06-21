BEGIN { 
  if( db != "") db=1;
  gsub(/\n/," ",libs);
  gsub(/\r/,"",libs);
  gsub(/\n/,"",libs);
  gsub(/[ ][ ]*/," ",libs);
  sub(/^ /,"",libs);
  sub(/ $/,"",libs);
  gsub(/([]])/," ] ",libs); 
  gsub(/([[])/," [ ",libs); 
  gsub(/([|])/," | ",libs); 
  if(db) print " LIBS='" libs "'";
  nl=split(libs,l," *"); 
  err= 0;
  sp=" ";
}
$2=="LD" { 
    if( ld[$1]=="" )  {
	a=$1;
	sub(/.* LD /,"");
	ld[a]=$0 sp; 
	if(db) print a,$0
	ex[a]=1;
    }}
$2=="INCLUDE" { if( inc[$1]=="" ) {
	a=$1;sub(/.* INCLUDE /,"");
	inc[a]=$0 sp;
	if(db)	print a,$0;
	ex[a]=1;
}}
END {
    msg=cpp sp;
    lvl=0;
    ok=1;
    skip=0; 
    k=0;
   
    for(i=1; i<=nl;++i)
    {
	lib=l[i];
 
	if(db) print " @@ '" lib "'", k,ok , skip , "err=" err, lvl,i,length(lib)
	if (length(lib)==0)
	{
	   if(db) print " empty "
	} 
	else if (lib=="[") 
	{ 
	    ncase=0;
	    skip=0;
	    lvl++;
	    ok=1;
	    k0=k;
	    lmis="";
	    if(lvl!=1)
		  lerr[err++]=" [ ... [ ... ";
	} 
	else if (lib=="]") 
	{ 
	    skip=ok ||skip;
	    if(skip==0) k=k0;
	    lvl--;
	    if( !skip && ncase) lerr[err++]=lmis  ;
	    ok=1;
	    skip=0;
	}
	else if(lib=="|") 
	{
	    ncase++;
	    if(db) print " |||| ",skip,ok,k,k0;
	    skip=ok || skip;
	    ok=!skip;
	    if(skip==0) k=k0;
	}
	else if( (skip==0) && (ok == 1))
	{
	    if(  ex[lib]!=1) {ok=0;
		if (lvl == 0) lerr[err++] = lib; 
		lmis=lmis " or " lib; }
	    else 
	    { 
		ln[++k]=lib;
		with[lib]=lvl;
	    }
	}
	if(db)  print "err=",err,"i=",i;
    }
    if(lvl !=0) 
	 lerr[err++]="no matching [|] ";

    if(db) 
    {
	printf "%s "," libs ::: " 
	for(j=1; j<=k; ++j)
	    printf "%s ",ln[j];
	print " err=" err,k,lvl;
    };
#  remove first  item  ..
    for (j=1; j<=k;++j)
	jln[ln[j]]=j;
  
    if(err==0) 
	for(j=1; j<=k; ++j)
	    if( jln[ln[j]]==j)
	    {
		lib = ln[j];
		if ( with[lib]) 
		    msg1="-DWITH_" lib sp ld[lib] sp inc[lib] sp;
		else
		    msg1= ld[lib] sp inc[lib] sp;
		
		msg =msg msg1; 
		if(db) print " ###" ,lib, err, first ,j, m, msg1;
	    }
    
    if(err) 
    {
	printf  "\t\t  MISSING lib "> "/dev/stderr" 
	for(i=0; i < err;++i)
	    printf "%s, ", lerr[i]> "/dev/stderr"; 
	print "         Check the WHERE-LIBRARYfiles " > "/dev/stderr";
	exit 1; 
    }
    else printf("%s\n",msg); 
}