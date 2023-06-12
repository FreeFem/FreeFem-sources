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
# SUMMARY : Makefile for FreeFEM
# LICENSE : LGPLv3
# ORG     : LJLL Universite Pierre et Marie Curie, Paris, FRANCE
# AUTHORS : ...
# E-MAIL  : ...

LIST_IDP=cobrameshcavity.idp cube.idp \
	DDM-funcs-v2.idp DDM-Schwarz-macro.idp \
	ExtractDofonBorder.idp ffddm_coarsemeshCS.idp \
	ffddm_functions.idp ffddm_geneoCS_3rdlvl.idp \
	ffddm_geneoCS.idp ffddm_parameters.idp \
	ffddm_partitioning.idp ffddm.idp \
	func_max.idp getARGV.idp \
	gsl.idp Heat3d.idp \
	macro_ddm_substructuring.idp macro_ddm.idp \
	MeshSurface.idp mortar-msh.idp \
	movemeshsmooth.idp MPIGMRESmacro.idp \
	MPIplot.idp buildmeshS.idp 

install-exec-local:
	$(mkinstalldirs) -m 755 "$(DESTDIR)$(ff_prefix_dir)/idp"
	$(INSTALL) -m 555 $(LIST_IDP) "$(DESTDIR)$(ff_prefix_dir)/idp"
