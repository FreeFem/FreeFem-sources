# ======================================================================
# Laboratoire Jacques-Louis Lions
# Université Pierre et Marie Curie-Paris6, UMR 7598, Paris, F-75005 France
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
# headeralh default=0 freefem make multipleauthors start=13/05/04 upmc

dist_pkgdata_DATA=$(DOCPDF)
documentation: $(DOCPDF) 
docps:  $(DOCPS) $(DOCPSGZ) 
docpdf: $(DOCPDF) 

# ALH - 16/10/13 - split GNU make directives into a separate makefile [[file:gnu.mak]] to suppress automake warning
# messages about non-portability

freefem++doc.pdf:
	$(MAKE) -f gnu.mak $@ EPSTOPDF="$(EPSTOPDF)"
freefem++doc.ps:
	$(MAKE) -f gnu.mak $@ EPSTOPDF="$(EPSTOPDF)"
freefem++doc.ps.gz:
	$(MAKE) -f gnu.mak $@ EPSTOPDF="$(EPSTOPDF)"
clean:
	$(MAKE) -f gnu.mak $@

# Local Variables:
# mode:makefile
# ispell-local-dictionary:"british"
# coding:utf-8
# End:
