# ======================================================================
# Modified by Antoine Le Hyaric
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
# headeralh default=0 freefem make modified start=19/03/10 upmc

include Makefile.inc

all-local: mumps

# Downloading and compiling mumps
# ------------------------------
DIRPKG=../pkg
SRCDIR=MUMPS_$(VERSION)
PACKAGE1=$(DIRPKG)/MUMPS_$(VERSION).tar.gz
INSTALL=../..
VERSION=4.10.0

mumps: FAIRE-$(VERSION)

$(SRCDIR)/FAIT: $(SRCDIR)/tag-tar
	cp Makefile.inc $(SRCDIR)
	cd $(SRCDIR);make d z
	touch $@
install: $(SRCDIR)/FAIT
	-mkdir -p ../include/libseq
	cp $(SRCDIR)/include/*.h ../include/
	cp $(SRCDIR)/libseq/*.h  ../include/libseq
	-mkdir -p ../lib
	cp $(SRCDIR)/lib/*.a ../lib/
	cp $(SRCDIR)/libseq/libmpiseqFREEFEM.a ../lib/

# ALH - WHERE was not built properly during a full compile from a
# clean source. This goal has been rewritten to solve this issue.

WHERE: 
	echo mumps LD -L@DIR@/lib   -ldmumpsFREEFEM -lzmumpsFREEFEM  -lmumps_commonFREEFEM -lpordFREEFEM \
		-lpthread  >../lib/WHERE.mumps
	echo mumps INCLUDE -I@DIR@/include  >> ../lib/WHERE.mumps

FAIRE-$(VERSION):
#
# ALH - the WHERE target must be built sequentially only after 'install'
#
	$(MAKE) install
	$(MAKE) WHERE 
	touch $@

Makefile.inc: ../../config.status	Makefile Makefile-mumps-$(VERSION).inc
	../../config.status  --file="Makefile.inc:Makefile-mumps-$(VERSION).inc"

$(SRCDIR)/$(INSTALL): $(SRCDIR)/tag-tar

$(SRCDIR)/tag-tar:$(PACKAGE1)
	tar xvzf $(PACKAGE1)
	patch -d MUMPS_$(VERSION)   -p 1    <MUMPS_$(VERSION).patch
	touch $(SRCDIR)/tag-tar

$(PACKAGE1):
	cd `dirname $@`; $(WGET)   http://graal.ens-lyon.fr/MUMPS/`basename $(PACKAGE1)`

clean-local:
	-cd $(SRCDIR) &&  $(MAKE) clean -C $(SRCDIR) 
	-rm Makefile.inc FAIRE* ../lib/WHERE.mumps
	-rm -rf ../include/*mumps*
	-rm -rf ../lib/lib*mumps* ../lib/libpord*.a ../lib/libmpiseq*.a
	-rm -rf $(SRCDIR)

clean: clean-local

# Local Variables:
# mode:makefile
# ispell-local-dictionary:"british"
# coding:utf-8
# End:
