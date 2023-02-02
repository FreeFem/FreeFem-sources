############################################################################
# This file is part of FreeFEM.                                            #
#                                                                          #
# FreeFEM is free software: you can redistribute it and/or modify          #
# it under the terms of the GNU Lesser General Public License as           #
# published by the Free Software Foundation, either version 3 of           #
# the License, or (at your option) any later version.                      #
#                                                                          #
# FreeFEM is distributed in the hope that it will be useful,               #
# but WITHOUT ANY WARRANTY; without even the implied warranty of           #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
# GNU Lesser General Public License for more details.                      #
#                                                                          #
# You should have received a copy of the GNU Lesser General Public License #
# along with FreeFEM. If not, see <http://www.gnu.org/licenses/>.          #
############################################################################
# SUMMARY : ...
# LICENSE : LGPLv3
# ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
# AUTHORS : ...
# E-MAIL  : ...

# To create statically linked executables (see configure.ac)
CXXLD = $(STATICTOOL) $(CXX)
bin_PROGRAMS = FreeFem++ FreeFem++-nw @FFGLUTPROG@
EXTRA_PROGRAMS = ffglut

# FFCS:visualization stream redirection
ffglut_SOURCES = ../femlib/Drawing.cpp \
  ../femlib/fem.cpp \
  ../fflib/ffapi.cpp \
  ../femlib/FQuadTree.cpp \
  ../femlib/GQuadTree.cpp \
  ../femlib/Mesh1dn.cpp \
  ../femlib/Mesh2dn.cpp \
  ../femlib/Mesh3dn.cpp \
  ../femlib/MeshSn.cpp \
  ../femlib/MeshLn.cpp \
  ../femlib/mshptg.cpp \
  ../Graphics/ffglut.cpp \
  ../Graphics/ffthreads.cpp \
  ../Graphics/ffthreads.hpp \
  ../Graphics/gggg.cpp

ffglut_DEPENDENCIES = ../libMesh/libMesh.a
ffglut_LDADD = ../libMesh/libMesh.a @LIBSGLUT@ @LIBSPTHREAD@
ffglut_LIBS =

FreeFem___nw_SOURCES = ../Graphics/sansrgraph.cpp ../mpi/parallelempi-empty.cpp ../fflib/ffapi.cpp ../fflib/compositeFESpace.cpp
FreeFem___SOURCES = ../Graphics/sansrgraph.cpp ../mpi/parallelempi-empty.cpp ../fflib/ffapi.cpp ../fflib/compositeFESpace.cpp
FreeFem___LDADD = ../lglib/liblg.a ../fflib/libff.a  @UMFPACKLIBS@ @ARPACKLIBS@ @BLASLIBS@ @FF_LIBS@ -lz
FreeFem___nw_LDADD = ../lglib/liblg.a ../fflib/libff.a  @UMFPACKLIBS@ @ARPACKLIBS@ @BLASLIBS@ @FF_LIBS@ -lz
FreeFem___nw_DEPENDENCIES = ../fflib/libff.a ../lglib/liblg.a 
FreeFem___DEPENDENCIES = ../fflib/libff.a ../lglib/liblg.a

#LDADD = ../lglib/liblg.a ../fflib/libff.a 

AM_CPPFLAGS = -I$(srcdir)/../fflib -I$(srcdir)/../Graphics -I$(srcdir)/../.. -I$(srcdir)/../femlib $(UMFPACK_CPPFLAGS) -DHAVE_ZLIB
