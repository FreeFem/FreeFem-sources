URL="http://ab-initio.mit.edu/nlopt/nlopt-2.2.4.tar.gz"
SRCDIR=nlopt-2.2.4
TARGZ=../pkg/$(SRCDIR).tar.gz
include Make.defs

all: FAIRE

FAIRE: ../pkg $(SRCDIR)/FAIT 

$(SRCDIR)/FAIT: $(TARGZ)
	tar zxvf $(TARGZ)
	cd $(SRCDIR); ./configure  --with-cxx --prefix="$(FFDOWNLOAD)" CXX="$(CXX)" CXXFLAGS="$(CXXFLAGS)" CC="$(CC)" CFLAGS="$(CFLAGS)" 
	cd $(SRCDIR); $(MAKE) install
	touch $(SRCDIR)/FAIT
	$(MAKE) WHERE
$(TARGZ):
	cd ../pkg;$(WGET) $(URL)
../pkg:
	mkdir $@
clean:
	rm -rf $(SRCDIR) Make.defs 
WHERE: 
	-if [ -f $(SRCDIR)/FAIT ] ; then \
	echo nlopt LD -L@DIR@/lib -lnlopt_cxx  >../lib/WHERE.nlopt ;\
	echo nlopt INCLUDE -I@DIR@/include >> ../lib/WHERE.nlopt ;\
	fi

Make.defs: ../../config.status	Makefile Make.inc
	../../config.status  --file="Make.defs:Make.inc"
