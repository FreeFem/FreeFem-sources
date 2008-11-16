# Makefile using Automake + Autoconf
# ----------------------------------
# $Id$
CXXLD=$(STATICTOOL) $(CXX)

bin_PROGRAMS=$(MEDITPROG)
EXTRA_PROGRAMS=ffmedit
ffmedit_SOURCES=  \
animat.c		cube.c			image.c			listnum.c		mouse.c			persp.c			stream.c		zaldy2.c \
bbfile.c		dlists.c		inmsh2.c		material.c		normal.c		picking.c		tiles.c \
camera.c		ellipse.c		inout.c			medit.c			param.c			psfile.c		transform.c \
cenrad.c		geometry.c		inout_morice.c		menus.c			parsar.c		scene.c			util.c \
clip.c			gisfil.c		inout_popenbinaire.c	mesh.c			parsop.c		scissor.c		vector.c \
clipvol.c		hash.c			items.c			mlists.c		particle.c		sftcpy.c		view.c \
critip.c		ilists.c		keyboard.c		morphing.c		path.c			status.c		zaldy1.c

LDADD= ../libMesh/libMesh.a $(LIBSGLUT)
INCLUDES=-I../libMesh
BUILT_SOURCES=compil.date
compil.date: 
	echo "#define COMPIL " '"' `date` '(with ff++ $(VERSION))''"' > compil.date
.PHONY: compil.date

