# egrep "^int|^double" /opt/local/include/gsl/gsl_sf* 
# egrep "^int|^double|^long" /usr/local/include/gsl/gsl_{cdf,randist,sf_*}.h | egrep -v ',$'  
#  export  LANG=C
BEGIN {FS="[(),]";
    FLOAD="gsl"
    LANG="C";
    edp=FLOAD ".idp";
    tex=FLOAD ".tex";
    print "\\begin{itemize}\n" > tex
    c="\"";
    pp=c "(" c;
    print " /*  ";
    print "// file create:  awk -f gsl.awk  gsl_list  > ff_gsl_awk.hpp " >edp
    print "load",c FLOAD c  > edp
    print "gslrng ffrng;" >edp
    print " gslabortonerror=0; " >edp
    #  convertion de type gsl .. ff++ 
    
    Tff["double"]="double";
    Tff["const double"]="double";
    Tff["int"]="long";
    Tff["const size_t"]="long";
    Tff["unsigned int"]="long";
    Tff["const int"]="long";
    Tff["const unsigned int"]="long";
    Tff["const double"]="double";
    Tff["const gsl_rng *"] = "gsl_rng **";
    Tff["gsl_mode_t"] = "long"; 
    Cff["const gsl_rng *"] = "*"; 
    V0["double"]= 0.55;
    V0["const double"]= 0.55; 
    V0["const gsl_rng *"] = "ffrng";
 }

function f0(f) { xx=f;gsub("_","",xx) ; return  xx  };
function f0tex(f) { xx=f;gsub("_","",xx) ; return  "\\texttt{" xx "}"  };
function ftex(f) { xx=f;gsub("_","\\_",xx) ; return  "\\texttt{" xx "}"  };
function ff(f) { return c f0(f) c };
function VV(t) { if ( V0[t] !=0) return V0[t]; else return 0;}
function ana1()
{
    ok = 0; 
    xx1=$1;
    nnn=split($1,xxx,":");
    if(nnn == 2) xx1= xxx[2]; 

    fonc=0; 
    split(xx1,fff,"[ \t]*");
    f=fff[2];
    R=fff[1]; 

    if( R=="int" || R =="double" ) ok=1;
    #print " xx: " xx1, R, f, ok ; 

}
function anatype(kkk)
{
    if( ok !=1) return 0; ; 
    v0=1;  
    t = "";
 
    nnn=split($(kkk+1),fff,"[ \t]*");
    i2 =1; 
    while( i2<nnn)
	{
	   if( t=="")
	    t = fff[i2];
	   else
	    t= t " " fff[i2]; 
	    i2++; 
	}
	tmp=match(fff[nnn],/.*\[\]/);
	if(tmp) t = t "[]"; 
	T[kkk]=t;
	if( Tff[t] == "") ok=0; 
	if( Tff[t] == "" && T00[t]=="") { print "//  -- missing type \"" t  "\" ";T00[t]=1;} 
    return ok; 
}

NF ==3 {
    ana1();
    #print NF, "******",ok
    anatype(1);
    if( ok == 1)
    {
      ok==2;
 	  lw="";
	  
	   lw = Tff[R] " " f "__(" Tff[T[1]] " const & x ) { return " f "( (" T[1] ")" Cff[T[1]] " x );}\n" ;
	  
	gg = "\n   Global.Add(" ff(f) "," pp ",new OneOperator1_<" Tff[R] "," Tff[T[1]] ">( " f "__)); "
	cw = cw  lw;
	cm = cm  gg;
	v0 = VV(T[1]);
	print "cout << " c  f "("v0") =  " c " << " f0(f) "("v0")  << endl; "> edp
	print "\\item " ftex(f) "(a) $\\quad \\mapsto\\quad $ "  f0tex(f) "(a) " > tex; 

	}
}

NF ==4 {
    ana1();
    anatype(1);
    anatype(2);
    if( ok == 1)
    {
      ok==2;
 	  lw="";
	  
	   lw = Tff[R] " " f "__(" Tff[T[1]] " const & x , " Tff[T[2]] " const & y )" \
	   "{ return " f "( (" T[1] ")" Cff[T[1]] " x , (" T[2] ")" Cff[T[2]] " y );}\n" ;
	  
	gg = "\n   Global.Add(" ff(f) "," pp ",new OneOperator2_<" Tff[R] "," Tff[T[1]] "," Tff[T[2]] ">( " f "__)); "
	cw = cw  lw;
	cm = cm  gg;
	v0 = VV(T[1]);
	v1 = VV(T[2]);
	
	print "cout << " c  f "(" v0 ", "v1 ") =  " c " << " f0(f) "(" v0 ",",v1 ")  << endl; "> edp
	print "\\item " ftex(f) "(a,b) $\\quad \\mapsto\\quad $ "  f0tex(f) "(a,b) " > tex; 
	}
}

NF ==5 {
    ana1();
    anatype(1);
    anatype(2);
    anatype(3);
    if( ok == 1)
    {
      ok==2;
 	  lw="";
	  
	   lw = Tff[R] " " f "__(" Tff[T[1]] " const & x , "  Tff[T[2]] " const & y , " Tff[T[3]] " const & z )" \
	   "{ return " f "( (" T[1] ")" Cff[T[1]] " x , (" T[2] ")" Cff[T[2]] " y , (" T[3] ")" Cff[T[3]] " z );}\n" ;
	  
	gg = "\n   Global.Add(" ff(f) "," pp ",new OneOperator3_<" Tff[R] "," Tff[T[1]] "," Tff[T[2]] "," Tff[T[3]] \
	   ">( " f "__)); ";
	cw = cw  lw;
	cm = cm  gg;
	v0 = VV(T[1]);
	v1 = VV(T[2]);
	v2 = VV(T[3]);
	
	print "cout << " c  f "(" v0 "," v1 "," v2 ") =  " c " << " f0(f) "(" v0 "," v1 "," v2 ")  << endl; "> edp
	print "\\item " ftex(f) "(a,b,c) $\\quad \\mapsto\\quad $ "  f0tex(f) "(a,b,c) " > tex; 

	}
}


ok !=1 { 
   print " missing:",NF," ", f , " -> " $0;
}
END {   print " */ "
	print "/*****************/";
	print "/*****************/";
    print "\\end{itemize}\n" > tex
	print cw ;

	print "/*****************/";
	print "/*****************/";
	print " void init_gsl_sf() { \n"
 	print cm ; 
	print " } ";

	print "/*****************/";
	print "/*****************/";

    }