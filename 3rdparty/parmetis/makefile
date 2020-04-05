# Downloading and compiling extra libraries
# -----------------------------------------

include Makefile.in

all-local: parmetis

# Downloading and compiling Tetgen
# ------------------------------
# http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-4.0.tar.gz
# Parmetis information
DIRPKG=../pkg
SRCDIR=parmetis-$(parmetis_VERSION)
PACKAGE=$(DIRPKG)/ParMetis-$(parmetis_VERSION).tar.gz
SERVER=http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis

INSTALL=../..
parmetis_VERSION=4.0.3

parmetis: FAIRE
	test -z "$(MPIRUN)" ||$(MAKE) FAIRE-mpi
$(SRCDIR)/FAIT:
	$(MAKE) install.done
	test -z "$(MPIRUN)" ||$(MAKE) $(SRCDIR)/FAIT-mpi 
	-test -z "$(MPIRUN)" &&  rm $(SRCDIR)/FAIT-mpi
	touch $(SRCDIR)/FAIT
$(SRCDIR)/FAIT-mpi:
	make install-mpi.done
	touch $(SRCDIR)/FAIT-mpi

install.done:$(SRCDIR)/tag-compile
	-mkdir -p ../include
	-mkdir -p ../lib
	#test -z "$(MPIRUN)" || cp $(SRCDIR)/libparmetis.a ../lib
	cp $(SRCDIR)/libmetis.a ../lib
	cp $(SRCDIR)/include/*.h ../include
	cp $(SRCDIR)/metis/include/*.h ../include
	touch $@
install-mpi.done:$(SRCDIR)/tag-compile-mpi
	cp $(SRCDIR)/libparmetis.a ../lib
	touch $@

FAIRE: $(SRCDIR)/FAIT 
	$(MAKE) WHERE
	test -z "$(MPIRUN)" || $(MAKE) WHERE-mpi
	touch FAIRE
FAIRE-mpi: $(SRCDIR)/FAIT-mpi 
	$(MAKE) WHERE-mpi
	touch FAIRE-mpi

METISCFLAGS=$(FLAGS_REGEX_H) $(CFLAGS)  -I../GKlib -I../include  -I. -D__thread= 
PMETISCFLAGS=$(MPI_INCLUDE) $(FLAGS_REGEX_H) $(CFLAGS)  -I../metis/GKlib -I../metis/include   -I../include -I. -D__thread= 
$(SRCDIR)/tag-compile: $(SRCDIR)/tags
	make -C $(SRCDIR)/metis/libmetis $(OBJS_METIS) "CC=$(CC)" "CFLAGS=$(METISCFLAGS)"
	make -C $(SRCDIR)/metis/GKlib $(OBJS_GK) "CC=$(CC)" CFLAGS="$(METISCFLAGS) -Dmetis_EXPORTS"
#	test -z "$(MPIRUN)" ||  make -C $(SRCDIR)/libparmetis $(OBJS_PARMETIS) "CC=$(CC)" "CFLAGS=$(PMETISCFLAGS)"
	-mkdir $(SRCDIR)/lib
	cd $(SRCDIR); $(AR)  libmetis.a  metis/libmetis/*.o metis/GKlib/*.o
#	test -z "$(MPIRUN)" ||(cd $(SRCDIR); $(AR)  libparmetis.a  libparmetis/*.o )
	touch $(SRCDIR)/tag-compile
$(SRCDIR)/tag-compile-mpi: $(SRCDIR)/tags
	 make -C $(SRCDIR)/libparmetis $(OBJS_PARMETIS) "CC=$(CC)" "CFLAGS=$(PMETISCFLAGS)"
	-mkdir $(SRCDIR)/lib
#	cd $(SRCDIR); $(AR)  libmetis.a  metis/libmetis/*.o metis/GKlib/*.o
	(cd $(SRCDIR); $(AR)  libparmetis.a  libparmetis/*.o )
	touch $(SRCDIR)/tag-compile-mpi



WHERE: install.done
	-if [ -f $(SRCDIR)/FAIT ] ; then \
		echo metis LD -L@DIR@/lib  -lmetis  >$(SRCDIR)/$(INSTALL)/lib/WHERE.metis ;\
		echo metis INCLUDE -I@DIR@/include   >> $(SRCDIR)/$(INSTALL)/lib/WHERE.metis ;\
	fi
WHERE-mpi: install-mpi.done
	-if [ -f $(SRCDIR)/FAIT-mpi ] ; then \
		echo parmetis LD -L@DIR@/lib -lparmetis -lmetis  >$(SRCDIR)/$(INSTALL)/lib/WHERE.parmetis ;\
		echo parmetis INCLUDE -I@DIR@/include   >> $(SRCDIR)/$(INSTALL)/lib/WHERE.parmetis ;\
	fi


OBJS_METIS=auxapi.o	coarsen.o	fm.o		graph.o	kwayrefine.o	minconn.o	options.o	separator.o	timing.o \
balance.o	compress.o	fortran.o	initpart.o	mcutil.o	mincover.o	parmetis.o	sfm.o		util.o \
bucketsort.o	contig.o	frename.o	kmetis.o	mesh.o	mmd.o		pmetis.o	srefine.o	wspace.o \
checkgraph.o	debug.o	gklib.o	kwayfm.o	meshpart.o	ometis.o	refine.o	stat.o

OBJS_GK=b64.o		error.o	fs.o		graph.o	itemsets.o	omp.o		random.o	sort.o	tokenizer.o \
blas.o	evaluate.o	getopt.o	htable.o	mcore.o	pdb.o		rw.o		string.o	util.o \
csr.o		fkvkselect.o	gkregex.o	io.o		memory.o	pqueue.o	seq.o		timers.o 

OBJS_PARMETIS=akwayfm.o		diffutil.o		initpart.o		move.o		renumber.o		wave.o \
ametis.o		frename.o		kmetis.o		msetup.o		rmetis.o		weird.o  \
balancemylink.o	gkmetis.o		kwayrefine.o		node_refine.o		selectq.o		wspace.o \
comm.o		gkmpi.o		match.o		ometis.o		serial.o		xyzpart.o \
csrmatch.o		graph.o		mdiffusion.o		pspases.o		stat.o \
ctrl.o		initbalance.o		mesh.o		redomylink.o		timer.o \
debug.o		initmsection.o	mmetis.o		remap.o		util.o 



# FFCS: patch is necessary for metis 4.0, but not for 4.0.3
$(SRCDIR)/tags: $(PACKAGE)
	tar xzf $(PACKAGE)
	patch -p1  <parmetis-4.0.3.patch
	touch $(SRCDIR)/tags




compile.done-old:$(SRCDIR)/tags Makefile.in
	-mkdir -p "../lib/parmetis"
	-mkdir -p "../include/parmetis"
	cp Makefile.in $(SRCDIR)
# 	FFCS - we need to know when errors occur
	cd $(SRCDIR)/METISLib && $(MAKE) 
	cd $(SRCDIR)/ParMETISLib && $(MAKE) 
	touch $@

# FFCS - simpler makefile
install.done-old:compile.done
#	cp $(SRCDIR)/libparmetis.a $(SRCDIR)/$(INSTALL)/lib
#	cp $(SRCDIR)/libmetis.a $(SRCDIR)/$(INSTALL)/lib
	cp $(SRCDIR)/parmetis.h ../include/parmetis
	cp $(SRCDIR)/METISLib/*.h ../include/parmetis
	sed  's;../parmetis.h;parmetis.h;' <$(SRCDIR)/METISLib/metis.h >../include/parmetis/metis.h
#
#	FFCS - 23/5/12 - cannot keep name libmetis.a because it is identical to the library created by
#	[[file:../metis/Makefile]] and library path mechanisms at link time pick one for the other on MinGW.
	mv $(SRCDIR)/$(INSTALL)/lib/parmetis/libmetis.a $(SRCDIR)/$(INSTALL)/lib/parmetis/libmetis.a
#
	-cd $(SRCDIR)/Programs && make 
	touch $@

# FFCS - simpler makefile
WHERE.done:install.done
	echo parmetis LD -L@DIR@/lib -lparmetis -lmetis  >$(SRCDIR)/$(INSTALL)/lib/WHERE.parmetis ;
	echo parmetis INCLUDE -I@DIR@/include   >> $(SRCDIR)/$(INSTALL)/lib/WHERE.parmetis ;
	touch $@

Makefile.in: ../../config.status Makefile-parmetis.in
	../../config.status  --file="Makefile.in:Makefile-parmetis.in"



$(SRCDIR)/tags-old: $(PACKAGE)
	tar xzf $(PACKAGE)
	mv  $(SRCDIR)/Programs/Makefile $(SRCDIR)/Programs/Makefile-orig
	sed -e 's;BINDIR = .*$$;BINDIR = ../$(INSTALL)/bin;' \
            -e 's;../libparmetis.a;../$(INSTALL)/lib/parmetis/libparmetis.a;' \
            -e 's;../libmetis.a;../$(INSTALL)/lib/parmetis/libmetis.a;' \
            -e 's;-L[.][.];-L../$(INSTALL)/lib/parmetis;' \
          <$(SRCDIR)/Programs/Makefile-orig \
          >$(SRCDIR)/Programs/Makefile
	mv $(SRCDIR)/ParMETISLib/Makefile $(SRCDIR)/ParMETISLib/Makefile-orig
	sed   -e 's;../libparmetis.a;../$(INSTALL)/lib/parmetis/libparmetis.a;' \
	      -e 's;../libmetis.a;../$(INSTALL)/lib/parmetis/libmetis.a;' \
          <$(SRCDIR)/ParMETISLib/Makefile-orig \
          >$(SRCDIR)/ParMETISLib/Makefile
	mv $(SRCDIR)/METISLib/Makefile $(SRCDIR)/METISLib/Makefile-orig
	sed   -e 's;../libparmetis.a;../$(INSTALL)/lib/parmetis/libparmetis.a;' \
	      -e 's;../libmetis.a;../$(INSTALL)/lib/parmetis/libmetis.a;' \
          <$(SRCDIR)/METISLib/Makefile-orig \
          >$(SRCDIR)/METISLib/Makefile
	touch $(SRCDIR)/tags

$(PACKAGE):
	../getall -o ParMETIS -a

# FFCS - 23/5/12 - corrected bug in removing .a from ../lib
clean:
	-cd $(SRCDIR) && $(MAKE) realclean  -C  $(SRCDIR)
	-rm Makefile.in $(SRCDIR)/tags
	-rm -rf ../lib/parmetis
	-rm -rf ../include/parmetis
#	-rm -rf ../include/metis
	-rm lib/libmetis.* ib/libparmetis.*
	-rm include/metis.h include/parmetis.h 
	-rm -rf $(SRCDIR)
	-rm *.done
	-rm FAIRE* WHERE*
