# Downloading and compiling extra libraries
# -----------------------------------------

# $Id$

all-local: $(DOWNLOAD_ARPACK)
EXTRA_DIST=ARmake.m4

# Downloading and compiling ARPACK
# --------------------------------

ARPACKLIB=ARPACK/libarpack_ff++.a

arpack: $(ARPACKLIB)

$(ARPACKLIB): ARPACK 
	if [ -n "@FF_LAPACKdir@" ] ; then \
	$(F77) -c `echo $(FFLAGS)\ |sed -e s/-O.\*\ // ` ARPACK/LAPACK/dlamch.f -o ARPACK/LAPACK/dlamch.o; \
	fi; \
	cd ARPACK && make lib 

ARPACK: arpack96.tar.gz patch.tar.gz ARmake.m4
	gunzip -c arpack96.tar.gz | tar xvf -
	gunzip -c patch.tar.gz | tar xvf -
	m4 -DFF_BLASLIB="$(BLASLIB)" \
	   -DFF_FC="@F77@" \
           -DFF_FFLAGS="@FFLAGS@" \
	   -DFF_LAPACKdir="@FF_LAPACKdir@" \
           -DFF_LDFLAGS="@LDFLAGS@" \
           -DFF_HOME=`pwd`/ARPACK \
	    ARmake.m4 >ARPACK/ARmake.inc
	cd ARPACK/LAPACK;

arpack96.tar.gz:
	@WGET@ -N http://www.caam.rice.edu/software/ARPACK/SRC/arpack96.tar.gz

patch.tar.gz:
	@WGET@ -N http://www.caam.rice.edu/software/ARPACK/SRC/patch.tar.gz

clean-local:
	-rm -r arpack++.tar.gz arpack++
	-rm -r arpack96.tar.gz patch.tar.gz ARPACK
