# Downloading and compiling extra libraries
# -----------------------------------------

# $Id$

all-local: $(DOWNLOAD_BLAS)

# Downloading and compiling the Blas
# ----------------------------------

BLASLIB=ATLAS/lib/ff++/libatlas.a

blas: $(BLASLIB)

$(BLASLIB): atlas3.6.0.tar.bz2
	tar xjf atlas3.6.0.tar.bz2
	@echo In the following questions, \"Architecture name\" must be ff++
	@echo "(Press ENTER)"
	@read
	cd ATLAS && make
	cd ATLAS && make install arch=ff++

atlas3.6.0.tar.bz2:
	wget -N http://heanet.dl.sourceforge.net/sourceforge/math-atlas/atlas3.6.0.tar.bz2

clean-local:
	-rm atlas*.tar.bz2
	-rm -r ATLAS
