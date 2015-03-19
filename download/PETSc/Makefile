# ======================================================================
# Laboratoire Jacques-Louis Lions
# Université Pierre et Marie Curie-Paris6, UMR 7598, Paris, F-75005 France
# http://www.ljll.math.upmc.fr/lehyaric
# ======================================================================
# This file is part of Freefem++
# 
# Freefem++ is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of
# the License, or (at your option) any later version.
# 
# Freefem++ is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with Freefem++; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
# ======================================================================
# headeralh default=0 freefem make multipleauthors start=19/03/10 upmc

include Makefile.inc

all-local:

# Downloading and compiling mumps
# -------------------------------

DIRPKG=../pkg
VERSION=3.5.2
PACKAGE=../pkg/petsc-lite-$(VERSION).tar.gz 
SRCDIR=petsc-$(VERSION)


PETSc:
	@if [  -n "$(MPI_INCLUDE)" ]; then $(MAKE) WHERE ;\
	else echo " -- No  Petc , no MPI"; fi


WHERE:../lib/WHERE.PETSc

../lib/WHERE.PETSc:$(SRCDIR)/tag-make
	echo PETSc LD -L@DIR@/lib   -ldPETSc > $@
	echo PETSc INCLUDE -I@DIR@/include  >> $@



$(SRCDIR)/tag-make:$(SRCDIR)/tag-conf
	cd $(SRCDIR) && $(MAKE)
	touch $@

$(SRCDIR)/FAIT:$(SRCDIR)/tag-make 
	cd $(SRCDIR) && $(MAKE) install 
	touch $@


Makefile.inc:../../config.status Makefile Makefile-PETSc.inc
	../../config.status --file="Makefile.inc:Makefile-PETSc.inc"

$(SRCDIR)/tag-conf:$(SRCDIR)/tag-tar
	cd $(SRCDIR) &&  ./configure \
	--prefix=abs_top_builddir/download \
	 'CFLAGS=$(CFLAGS)' 'CC=$(CC)' 'FC=$(FC)' 'FCFLAGS=$(FCFLAGS)' \
        --with-mpi-include='$(MPI_INC_DIR)' --with-mpi-lib='$(MPI_LIB) $(MPI_LIBFC)' \
        --with-mpiexec='$(MPIRUN)' \
	--with-blas-lapack-lib='$(LAPACKLIBS) $(BLASLIBS)' \
	--download-scalapack --download-metis --download-ptscotch \
	--download-mumps --download-ml --download-hypre --download-parmetis \
	--download-superlu --download-pastix --download-suitesparse --download-fftw \
	touch $@

$(SRCDIR)/tag-tar:$(PACKAGE)
	tar xvzf $(PACKAGE)
	touch $@

clean:clean-local
clean-local:
	-cd $(SRCDIR) &&  $(MAKE) clean -C $(SRCDIR) 
	-rm Makefile.inc FAIRE* ../lib/WHERE.PETSc
	-rm -rf ../include/*PETSc*
	-rm -rf ../lib/lib*PETSc* 
	-rm -rf $(SRCDIR)
#	-rm $(PACKAGE1)
	-rm config.log *.done

# Local Variables:
# mode:makefile
# ispell-local-dictionary:"british"
# coding:utf-8
# End:
