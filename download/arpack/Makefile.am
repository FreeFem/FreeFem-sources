# Downloading and compiling extra libraries
# -----------------------------------------

# $Id$

all-local: $(DOWNLOAD_ARPACK)

# Downloading and compiling ARPACK
# --------------------------------

ARPACKLIB=ARPACK/libarpack_ff++.a

arpack: $(ARPACKLIB)

$(ARPACKLIB): arpack96.tar.gz patch.tar.gz
	gunzip -c arpack96.tar.gz | tar xvf -
	gunzip -c patch.tar.gz | tar xvf -
	cp ARmake.inc ARPACK
	cd ARPACK && make lib home=`pwd` BLASLIB="$(BLASLIB)" \
		FC="$(F77)" FFLAGS="$(FFLAGS)" LDFLAGS="$(LDFLAGS)"

arpack96.tar.gz:
	wget -N http://www.caam.rice.edu/software/ARPACK/SRC/arpack96.tar.gz

patch.tar.gz:
	wget -N http://www.caam.rice.edu/software/ARPACK/SRC/patch.tar.gz

clean-local:
	-rm -r arpack++.tar.gz arpack++
	-rm -r arpack96.tar.gz patch.tar.gz ARPACK
