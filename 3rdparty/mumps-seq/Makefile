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

# Downloading and compiling extra libraries
# -----------------------------------------

include Makefile.inc

all-local: mumps
mumps:../lib/WHERE.mumpsseq

# Downloading and compiling mumps
# -------------------------------

DIRPKG=../pkg
SRCDIR=MUMPS_$(VERSION)
PACKAGE1=$(DIRPKG)/MUMPS_$(VERSION).tar.gz
INSTALL=../..
VERSION=5.0.2

# FFCS - install and WHERE need to be done sequentially, even in parallel builds

../lib/WHERE.mumpsseq: install.done
	echo mumps_seq LD -L@DIR@/lib -ldmumpsFREEFEM-SEQ -lzmumpsFREEFEM-SEQ -lmumps_commonFREEFEM-SEQ	\
		-lpordFREEFEM-SEQ $(LIBSPTHREAD)  > ../lib/WHERE.mumpsseq
	echo mumps_seq INCLUDE -I@DIR@/include/mumps_seq >> ../lib/WHERE.mumpsseq
	echo libseq LD -L@DIR@/lib  -lmpiseqFREEFEM-SEQ >>../lib/WHERE.mumpsseq
	echo libseq INCLUDE -I@DIR@/include -I@DIR@/include/mumps_seq >> ../lib/WHERE.mumpsseq

install.done: $(SRCDIR)/FAIT
	-mkdir -p ../include/mumps_seq
	cp $(SRCDIR)/include/*.h ../include/mumps_seq
	cp $(SRCDIR)/libseq/*.h  ../include/mumps_seq
	-mkdir -p ../lib
	cp $(SRCDIR)/lib/*.a ../lib/
	cp $(SRCDIR)/libseq/libmpiseqFREEFEM-SEQ.a ../lib/
	touch $@

$(SRCDIR)/FAIT:Makefile.inc $(SRCDIR)/tag-tar
	cp Makefile.inc $(SRCDIR)
#
#	FFCS - 22/5/12 - Mumps has difficulties compiling d & z in parallel
#
	cd $(SRCDIR) && $(MAKE) d
	cd $(SRCDIR) && $(MAKE) z
	touch $@

Makefile.inc: ../../config.status Makefile Makefile-mumps-$(VERSION).inc
	../../config.status  --file="Makefile.inc:Makefile-mumps-$(VERSION).inc"

$(SRCDIR)/tag-tar:$(PACKAGE1)
	tar xzf $(PACKAGE1)
	touch $@

$(PACKAGE1):
	../getall -o MUMPS -a

clean:
	-rm Makefile.inc FAIRE* *~ 
	-rm -rf ../include/*mumps*
	-rm -rf ../lib/lib*mumps*SEQ* ../lib/libpord*SEQ.a ../lib/libmpiseq*SEQ.a ../lib/WHERE.mumpsseq
	-rm -rf $(SRCDIR)
	-rm *.done config.log

# Local Variables:
# mode:makefile
# ispell-local-dictionary:"british"
# coding:utf-8
# End:
