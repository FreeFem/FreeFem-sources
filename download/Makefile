# Downloading and compiling extra libraries
# -----------------------------------------

# $Id$

# Default is not to download anything (i.e. use local system-specific
# versions)

# Downloading and compiling the Blas
# ----------------------------------

BLASLIB=libxxx
BLASDIR=ATLAS
BLASPACK=atlas3.6.0.tar.bz2
BLASGET=http://heanet.dl.sourceforge.net/sourceforge/math-atlas/atlas3.6.0.tar.bz2

blas: $(BLASLIB)
$(BLASLIB): $(BLASDIR)/Make.ff
	cd $(BLASDIR) && make install arch=ff++

$(BLASDIR)/Make.ff++: $(BLASDIR)
	cd $(BLASDIR) && make

$(BLASDIR): $(BLASPACK)
	tar xjf $(BLASPACK)

$(BLASPACK):
	wget $(BLASGET)

arpack:
	-mkdir arpack;
	cd arpack; \
	wget  http://www.caam.rice.edu/software/ARPACK/SRC/arpack96.tar.gz ;\
	wget  http://www.caam.rice.edu/software/ARPACK/SRC/patch.tar.gz  ;\
	wget  http://www.caam.rice.edu/software/ARPACK/ARPACK++/arpack++.tar.gz ; \
	gunzip -c arpack96.tar.gz | tar xvf - ;\
	gunzip -c patch.tar.gz | tar xvf - ;\
	gunzip -c arpack++.tar.gz| tar xvf - 
