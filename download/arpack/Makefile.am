# Downloading and compiling extra libraries
# -----------------------------------------

# $Id$

all-local: $(DOWNLOAD_ARPACK)
EXTRA_DIST=ARmake.m4 arpack-patch-lapack.tar.gz 

# Downloading and compiling ARPACK
# --------------------------------

# set in configure 
#ARPACKLIB=ARPACK/libarpack_ff++.a

arpack: $(ARPACKLIB)

$(ARPACKLIB): ARPACK/fait
	if [ -n '@FF_LAPACKdir@' ] ; then \
	$(F77) -c `echo $(FFLAGS)\ |sed -e s/-O.\*\ // ` ARPACK/LAPACK/dlamch.f -o ARPACK/LAPACK/dlamch.o; \
	fi; \
	cd ARPACK && make lib 
	$(AR) $(ARFLAGS) $(ARPACKLIB)  ARPACK/SRC/*.o  ARPACK/UTIL/*.o
	-if test -n "FF_LAPACKdir" ; then 	$(AR) $(ARFLAGS) $(LAPACKLIB)  ARPACK/LAPACK/*.o ; fi


ARPACK/fait: arpack96.tar.gz patch.tar.gz ARmake.m4 Makefile
	-rm -rf ARPACK
	gunzip -c arpack96.tar.gz | tar xf -
	gunzip -c patch.tar.gz | tar xf -
	gunzip -c arpack-patch-lapack.tar.gz | tar xf -
	for i in ARPACK/SRC/*.f ; do \
	    mv  $$i $$i.cpy; sed 's/call *second/call secnd2/' <$$i.cpy >$$i;rm $$i.cpy; done 
	for i in  ARPACK/UTIL/second.f;  do  \
	   mv  $$i $$i.cpy; cat $$i.cpy| sed 's/ SECOND *(/ secnd2(/'|grep -v EXTERNAL  >$$i;rm $$i.cpy; done 
	m4  -DFF_BLASLIB="$(BLASLIB)" \
	    -DFF_ARPACKLIB="$(ARPACKLIB)" \
            -DFF_LAPACKLIB="$(LAPACKLIB)" \
	    -DFF_FC="@F77@" \
            -DFF_FFLAGS="@FFLAGS@" \
	    -DFF_LAPACKdir='@FF_LAPACKdir@' \
            -DFF_LDFLAGS="@LDFLAGS@" \
            -DFF_HOME=`pwd`/ARPACK \
	    -DFF_SECOND="@FF_SECOND@" \
	    -DFF_AR="@AR@" \
	    -DFF_ARFLAGS="@ARFLAGS@" \
	    -DFF_RANLIB="@RANLIB@" \
	    ARmake.m4 >ARPACK/ARmake.inc
#
#	cp ARPACK/UTIL/Makefile ARPACK/Makefile_UTIL
#	sed  -e '1,$$s/second.o/$$(SECOND_O)/' <ARPACK/Makefile_UTIL  >ARPACK/UTIL/Makefile  
	touch ARPACK/fait
	cd ARPACK/LAPACK;

arpack96.tar.gz:
	@WGET@ -N http://www.caam.rice.edu/software/ARPACK/SRC/arpack96.tar.gz

patch.tar.gz:
	@WGET@ -N http://www.caam.rice.edu/software/ARPACK/SRC/patch.tar.gz

clean-local:
	-rm -r arpack++.tar.gz arpack++
	-rm -r arpack96.tar.gz patch.tar.gz ARPACK
