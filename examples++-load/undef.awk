# l=`otool -L MUMPS.dylib | awk ' NF>1 {printf("%s ", $1)}'`
# nm MUMPS.dylib | awk -f undef.awk -v l="$l" | c++filt -t
BEGIN { 
    if ( ! f) f="../src/mpi/FreeFem++-mpi"
    if( ! l ) l="MUMPS.dylib"; 
    FS="[ \t$]+"
    if(!ldd) lld="otool -L";
    lf = f " " l;
    {
	cldd=ldd " " lf;
	while (cldd " " f | getline >0)
	{
	    if(NF >1) ll=$1;
	    print " lib " ll; 
	    cmd="nm " ll;
	    # print cmd ;
	    while ( cmd | getline >0)
		if( $2 == "T" || $2 == "D" || $2 == "S" )   a[$3] =1; 
	    close(cmd) 
	}
	close(cldd)
    }
    {
	while ( cmd | getline >0)
	    if( $2 == "U"  && a[2] ==0 )  print $2; 
	close(cmd) 
    }
    exit;
}
