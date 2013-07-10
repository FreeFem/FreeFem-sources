# nm dfft.dylib | awk -f undef.awk | c++filt -t
BEGIN { 
    cmd="nm ../src/nw/FreeFem++ /usr/lib/libstdc++.6.dylib /usr/lib/libSystem.B.dylib /usr/local/lib/libgfortran.3.dylib /usr/lib/libc.dylib| sed 's/\$.*$//' ";
    while ( cmd | getline >0)
	if( $2 == "T" || $2 == "D" || $2 == "S" )   a[$3] =1; 
    close(cmd) 
 }

/ *U / { if( a[$2] ==0) print $2}