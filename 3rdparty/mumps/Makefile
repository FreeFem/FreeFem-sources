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

all-local: mumps

# Downloading and compiling mumps
# -------------------------------

DIRPKG=../pkg
SRCDIR=MUMPS_$(VERSION)
PACKAGE1=$(DIRPKG)/MUMPS_$(VERSION).tar.gz
VERSION=5.0.2

mumps:../lib/WHERE.mumps

$(SRCDIR)/FAIT: $(SRCDIR)/tag-tar
	cp Makefile.inc $(SRCDIR)
#
# ALH - 'd' and 'z' versions need to be made sequentially
#
	cd $(SRCDIR);make d
	cd $(SRCDIR);make z
	touch $@

# ALH - WHERE was not built properly during a full compile from a clean source. This goal has been rewritten to solve
# this issue.

../lib/WHERE.mumps:install.done
	echo mumps LD -L@DIR@/lib   -ldmumpsFREEFEM -lzmumpsFREEFEM  -lmumps_commonFREEFEM -lpordFREEFEM $(LIBSPTHREAD)  > $@
	echo mumps INCLUDE -I@DIR@/include  >> $@

install.done:$(SRCDIR)/FAIT
	-mkdir -p ../include/libseq
	cp $(SRCDIR)/include/*.h ../include/
	cp $(SRCDIR)/libseq/*.h  ../include/libseq
	-mkdir -p ../lib
	cp $(SRCDIR)/lib/*.a ../lib/
	cp $(SRCDIR)/libseq/libmpiseqFREEFEM.a ../lib/
	touch $@

FAIRE-$(VERSION):WHERE
	touch FAIRE-$(VERSION)

# ALH - "make d" and "make z" need to be done sequentially
$(SRCDIR)/FAIT:$(SRCDIR)/tag-tar Makefile.inc
	cp Makefile.inc $(SRCDIR)
	cd $(SRCDIR) && $(MAKE) d
	cd $(SRCDIR) && $(MAKE) z
	touch $@

Makefile.inc:../../config.status Makefile Makefile-mumps-$(VERSION).inc
	../../config.status --file="Makefile.inc:Makefile-mumps-$(VERSION).inc"

$(SRCDIR)/tag-tar:$(PACKAGE1)
	tar xzf $(PACKAGE1)
	#patch -d MUMPS_$(VERSION) -p 1 < MUMPS_$(VERSION).patch
	touch $(SRCDIR)/tag-tar

clean:clean-local
clean-local:
#FH	-cd $(SRCDIR) &&  $(MAKE) clean -C $(SRCDIR) 
	-rm Makefile.inc FAIRE* ../lib/WHERE.mumps
	-rm -rf ../include/*mumps*
	-rm -rf ../lib/lib*mumpsFREEFEM* ../lib/libpordFREEFEM.a ../lib/libmpiseqFREEFEM.a 
	-rm -rf $(SRCDIR)
#	-rm $(PACKAGE1)
	-rm config.log *.done

# Local Variables:
# mode:makefile
# ispell-local-dictionary:"british"
# coding:utf-8
# End:
