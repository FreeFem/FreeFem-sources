# Downloading and compiling extra libraries
# -----------------------------------------

# $Id$

SRCLIB=ATLAS/lib/ff++/libatlas.a
BINLIB=WinNT_PPRO256/libatlas.a

all-local: $(DOWNLOAD_BLAS)
source: $(SRCLIB)
binary: $(BINLIB)

# Downloading and compiling the Blas
# ----------------------------------

$(SRCLIB): atlas3.6.0.tar.bz2
	tar -x --bzip -f atlas3.6.0.tar.bz2
	@echo In the following questions, \"Architecture name\" must be ff++
	@echo "(Press ENTER)"
	@read
	cd ATLAS && make
	cd ATLAS && make install arch=ff++

atlas3.6.0.tar.bz2:
	wget -N http://heanet.dl.sourceforge.net/sourceforge/math-atlas/atlas3.6.0.tar.bz2

$(BINLIB): atlas321_WinNT_PPRO.zip
	unzip atlas321_WinNT_PPRO.zip
	touch $(BINLIB)

atlas321_WinNT_PPRO.zip:
	wget -N http://www.netlib.org/atlas/archives/windows/atlas321_WinNT_PPRO.zip

clean-local:
	-rm atlas*.tar.bz2
	-rm -r ATLAS
	-rm atlas321_WinNT_PPRO.zip
	-rm -r WinNT_PPRO256
