# l=`otool -L MUMPS.dylib | awk ' NF>1 {printf("%s ", $1)}'`
#  awk -f undef.awk -v l MUMPS.dylib -v f=../src/mpi/FreeFem++-mpi
BEGIN { 
    so=".dylib";
    if (  f==0) f="../src/mpi/FreeFem++-mpi"
    if( l==0 ) l="MUMPS.dylib"; 
    FS="[ \t$]+"
    if(ldd==0) ldd="otool -L ";
    print " cmmd : ",lld ", " l ", "f "." 
    nll=1;
    noo=O;
    lln[0]=l;
    lll[l]=0;
    while(noo < nll)
    {
	cldd=ldd " " lln[noo++]  " " f;
	#print " cldd = " cldd; 
	add=0; 
	while (cldd | getline >0)
	{
	   # print $2" ..." , $0;
	    ll=$2;
	    if( match(ll,"[.]dylib$")>0 && lll[ll] == 0 )
	    {
		print nll, ll;
		lln[nll++]=ll;
		lll[ll]=nll;
	    }
	}
	close(cldd)
    }
    print nll " " 
    for( il = 1; il < nll; il++)
    {
	ll=lln[il];
	cmd="nm " ll;
	## print " .... "cmd ;
	while ( cmd | getline >0)
	{
	    if( $2 == "T" || "t" == $2 || $2 == "D" || $2 == "S" )   { if(d) print " def :" $3 "."; a[$3]++;} 
	}
	close(cmd) 
    }
    
    {
	cmd="nm " l;
	while ( cmd | getline >0)
	    if( $2 == "U"  && a[$3] ==0 )  print $0, " Undef ???  " ;
	close(cmd) 
    }
    exit;
}
