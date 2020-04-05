URL="http://ab-initio.mit.edu/nlopt/nlopt-2.2.4.tar.gz"
SRCDIR=nlopt-2.2.4
TARGZ=../pkg/$(SRCDIR).tar.gz
-include Make.defs

all: FAIRE

FAIRE: ../pkg $(SRCDIR)/FAIT 

$(SRCDIR)/FAIT: $(TARGZ)
	tar zxf $(TARGZ)
#
#	FFCS - 30/11/11 - "--disable-dependency-tracking": dependencies pose problem on mixed Cygwin/MinGW setups
#	because of slashes and backslashes
#
# 	FFCS - 22/10/12 - Fred - "--without-octave": patch for one machine that has octave installed on it: building of
#	nlopt fails because it wants to install the octave components.  Not a bad idea in general but probably not when
#	building ffcs.
#
	#remove c++11 flags for windfows  FH ..
	cxxflags=`echo $(CXXFLAGS) |sed -e 's/[-]std=c[+][+]11//g'` ; \
	cd $(SRCDIR) && ./configure --disable-dependency-tracking --with-cxx --prefix="$(FFDOWNLOAD)" CXX="$(CXX)"	\
	CXXFLAGS="$(CXXFLAGS)" CC="$(CC)" CFLAGS="$(CFLAGS)" CPP="$$cxxflag" "--without-threadlocal" "--without-octave"
#
#	FFCS - 30/11/11 - we need to know when the building process breaks
#
	cd $(SRCDIR) && $(MAKE) install
	touch $(SRCDIR)/FAIT
	$(MAKE) WHERE
../pkg:
	mkdir $@

$(TARGZ):
	../getall -o NLopt -a

clean-local:
	rm -rf $(SRCDIR) Make.defs 
	rm -f ../lib/libnlopt_cxx*  ../lib/WHERE.nlopt 
	rm -f ../include/nlopt*

clean:clean-local

WHERE: 
	-if [ -f $(SRCDIR)/FAIT ] ; then \
	echo nlopt LD -L@DIR@/lib -lnlopt_cxx  >../lib/WHERE.nlopt ;\
	echo nlopt INCLUDE -I@DIR@/include >> ../lib/WHERE.nlopt ;\
	fi

Make.defs: ../../config.status	Makefile Make.inc
	../../config.status  --file="Make.defs:Make.inc"
