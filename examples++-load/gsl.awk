# egrep "^int|^double" /opt/local/include/gsl/gsl_sf* 
BEGIN {FS="[(),]";
    LANG="C";
    edp="gsl1.edp";
    c="\"";
    pp=c "(" c;
    print " /*  ";
    print "load",c "gsl" c  > edp
}

function f0(f) { xx=f;gsub("_","",xx) ; return  xx  };
function ff(f) { return c f0(f) c };
NF ==3 {
    split($1,fff,"[ \t]*");
    f=fff[2];
    split($2,fff,"[ \t]*");
    tt=0;
    i2=1
    ok=1
    nn=0;
    v0=1;
	t="";
	while(ok==1 && nn++<5)
	{
	    # print " ---- " , fff[i2] , i2 
	    if(fff[i2] == "unsigned" ) { t= t " " fff[i2] ;}
	    else if(fff[i2] == "int" ) { t= t " " fff[i2] ;tt="long";}
	    else if(fff[i2] == "double" ) { t= t " " fff[i2];tt="double"; v0=0.55 ;}
	    else if(fff[i2] == "const" ) { t= t " " fff[i2] ;}
	    else {i2--;ok=0;}
	    i2++;
	}
	p1=fff[i2];
	#print f, "  " , t , " " , tt, " " , p1;
	lw="";
	if(tt!=0){
	if(tt=="long")
	{
	    lw ="double " f "( long x ) { return " f "( (" t ") x );}\n" ;
	}
	gg = "\n   Global.Add("ff(f) "," pp",new OneOperator1<double ," tt ">( " f " )); "
	cw = cw  lw;
	cm = cm  gg;
	print "cout << " c  f "("v0") =  " c " << " f0(f) "("v0")  << endl; "> edp
	}
	else
	    print " missng "f,t,tt
}
NF > 4 { print " minssing " $0 ; }
NF ==4 {
    split($1,fff,"[ \t]*");
    f=fff[2];
    
    split($2,fff,"[ \t]*");
    tt=0;
    i2=1
    ok=1
    nn=0;
    t="";
    while(ok==1 && nn++<5)
    {
	# print " ---- " , fff[i2] , i2 
	if(fff[i2] == "unsigned" ) { t= t " " fff[i2] ;}
	else if(fff[i2] == "int" ) { t= t " " fff[i2] ;tt="long";}
	else if(fff[i2] == "double" ) { t= t " " fff[i2];tt="double" ;}
	else if(fff[i2] == "const" ) { t= t " " fff[i2] ;}
	else {i2--;ok=0;}
	i2++;
    }
    
    p1=fff[i2];
    tt1=tt;
    t1=t;
    
    split($3,fff,"[ \t]*");
    tt=0;
    i2=1
    ok=1
    nn=0;
    t="";
    if (fff[i2]==" ") i2++;
    if (fff[i2]=="") i2++;
    if (fff[i2]=="  ") i2++;
    while(ok==1 && nn++<5)
    {
	# print " ---- " , fff[i2] , i2 
	if(fff[i2] == "unsigned" ) { t= t " " fff[i2] ;}
	else if(fff[i2] == "int" ) { t= t " " fff[i2] ;tt="long";}
	else if(fff[i2] == "gsl_mode_t" ) { t= t " " fff[i2] ;tt="long";}
	else if(fff[i2] == "double" ) { t= t " " fff[i2];tt="double" ;}
	else if(fff[i2] == "const" ) { t= t " " fff[i2] ;}
	else {i2--;ok=0;}
	i2++;
    }
    
    p2=fff[i2];

#    print f, "  " , t , " " , tt, " " , p1, p3;
    lw="";
    if(tt && tt1 ) 
    {
	if(tt1=="long" || tt="long")  
	{
	    lw ="double " f "( " tt1  " x, " tt " y   ) { return " f "( (" t1 ") x, (" t ")  y  );}\n" ;
	}
	gg = "\n   Global.Add(" ff(f) "," pp ",new OneOperator2<double," tt1 "," tt " >( " f " )); ";
	cw = cw  lw;
	cm = cm  gg;
    }
    else
	print NF, "  ", tt1,",",tt," --" t1," ",t " ::: " " missing -- ", f , " (" , $2 "," $3  ") " fff[1] "," fff[2],",", fff[3] 
}

END {   print " */ "
	print "/*****************/";
	print "/*****************/";

	print cw ;

	print "/*****************/";
	print "/*****************/";
	print " void init_gsl_sf() { \n"
 	print cm ; 
	print " } ";

	print "/*****************/";
	print "/*****************/";

    }