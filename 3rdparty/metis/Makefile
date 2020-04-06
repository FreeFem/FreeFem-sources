# Downloading and compiling extra libraries
# -----------------------------------------

include Makefile.in

all-local: metis

# Downloading and compiling Tetgen
# ------------------------------
# http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-4.0.tar.gz
# Metis information
DIRPKG=../pkg
SRCDIR=metis-$(metis_VERSION)
PACKAGE=$(DIRPKG)/metis-$(metis_VERSION).tar.gz
SERVER=http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis
#//http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/OLD
INSTALL=../..

# FFCS - 14/11/11 - version 4.0.3 is not available from netlib anymore
metis_VERSION=5.1.0

metis: FAIRE


$(SRCDIR)/FAIT:
	$(MAKE) install
	touch $(SRCDIR)/FAIT

install:$(SRCDIR)/tag-compile
#	cd $(SRCDIR)/Programs ;$(MAKE)
#	cd $(SRCDIR) ;  $(MAKE)  -C  install
	-mkdir -p ../include
	-mkdir -p ../lib
	cp $(SRCDIR)/libmetis.a ../lib
	cp $(SRCDIR)/include/*.h ../include


FAIRE: $(SRCDIR)/FAIT
	$(MAKE) WHERE
	touch FAIRE

Makefile.in: ../../config.status Makefile-metis.in
	../../config.status  --file="Makefile.in:Makefile-metis.in"

$(SRCDIR)/tag-compile: $(SRCDIR)/tags
	#cd $(SRCDIR)/libmetis ;
	make -C $(SRCDIR)/libmetis $(OBJS) 'CC=$(CC)' CFLAGS='$(CFLAGS) -I../GKlib -I../include  -I.'
	make -C $(SRCDIR)/GKlib $(OBJS_GK) 'CC=$(CC)' CFLAGS='$(CFLAGS) -I../GKlib -I../include  -I. -Dmetis_EXPORTS -D_GNU_SOURCE'
	mkdir $(SRCDIR)/lib
	cd $(SRCDIR); $(AR)  libmetis.a  libmetis/*.o GKlib/*.o
	touch $(SRCDIR)/tag-compile



WHERE:
	-if [ -f $(SRCDIR)/FAIT ] ; then \
	echo metis LD -L@DIR@/lib -lmetis  >$(SRCDIR)/$(INSTALL)/lib/WHERE.metis ;\
	echo metis INCLUDE -I@DIR@/include >> $(SRCDIR)/$(INSTALL)/lib/WHERE.metis ;\
	fi


OBJS=auxapi.o	coarsen.o	fm.o		graph.o	kwayrefine.o	minconn.o	options.o	separator.o	timing.o \
balance.o	compress.o	fortran.o	initpart.o	mcutil.o	mincover.o	parmetis.o	sfm.o		util.o \
bucketsort.o	contig.o	frename.o	kmetis.o	mesh.o	mmd.o		pmetis.o	srefine.o	wspace.o \
checkgraph.o	debug.o	gklib.o	kwayfm.o	meshpart.o	ometis.o	refine.o	stat.o

OBJS_GK=b64.o		error.o	fs.o		graph.o	itemsets.o	omp.o		random.o	sort.o	tokenizer.o \
blas.o	evaluate.o	getopt.o	htable.o	mcore.o	pdb.o		rw.o		string.o	util.o \
csr.o		fkvkselect.o	gkregex.o	io.o		memory.o	pqueue.o	seq.o		timers.o



# FFCS: patch is necessary for metis 4.0, but not for 4.0.3
$(SRCDIR)/tags: $(PACKAGE)
	tar xzf $(PACKAGE)
	patch -p0  <metis-5.1.patch
	touch $(SRCDIR)/tags

$(PACKAGE):
	../getall -o METIS -a

clean-local:
	-cd $(SRCDIR) && $(MAKE) realclean  -C $(SRCDIR)
	-rm -rf config.log ../lib/metis ../lib/WHERE.metis
# remove metis if no parmetis 
	-test -f ../lib/libparmetis.a || rm ../lib/libmetis.a 
	-test -f ../include/parmetis.h || rm ../include/metis.h
	-rm -rf ../include/metis
	-rm -rf $(SRCDIR)
	-rm FAIRE FAIT

# FFCS - 23/5/12 - $(SRCDIR)/$(INSTALL) is meaningless if $(SRCDIR) does not exist
clean: clean-local
