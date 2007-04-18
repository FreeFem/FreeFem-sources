# Downloading and compiling extra libraries
# -----------------------------------------

# $Id$

all-local: $(DOWNLOAD_FLTK)

# Downloading and compiling FLTK
# ------------------------------

# FLTK information
SRCDIR=fltk-$(FLTK_VERSION)
PACKAGE=fltk-$(FLTK_VERSION)-source.tar.bz2
SERVER=http://switch.dl.sourceforge.net/sourceforge/fltk
INSTALL=../..

fltk: $(SRCDIR)/FF

# Using system JPEG and PNG lib is a problem under Cygwin
$(SRCDIR)/FFF: $(SRCDIR)
	cd $(SRCDIR) && ./configure --enable-threads --enable-localjpeg	\
		--enable-localpng $(FLTK_CONFIG_PARAM) --disable-gl	\
		@FLTK_XONLY@ --prefix=`pwd`/$(INSTALL) CXX=$(CXX) CC=$(CC)
	touch $(SRCDIR)/FFF
	cd $(SRCDIR); \
	if grep '#define'  config.h|egrep -q "U32|U64" ; then  \
	 echo " fltk conf OK ??? FH ";  \
	else \
	 echo " fltk configure  bug (on windows) ????  add #define U32 unsigned long (bofbof) " ;\
	 echo "#define U32 unsigned long" >> config.h ; \
	 echo "stop the make process ; just redo make to continue " ; \
	 exit 1; \
	 fi

$(SRCDIR)/FF: $(SRCDIR)/FFF
	cd $(SRCDIR) && make
	cd $(SRCDIR) && make install
	touch FF
$(SRCDIR): $(PACKAGE)
	tar xvjf $(PACKAGE)
	touch $(SRCDIR)
$(PACKAGE):
	@WGET@ -N $(SERVER)/$(PACKAGE)

clean-local:
	-rm -r fltk-*
