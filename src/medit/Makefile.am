# Makefile using Automake + Autoconf
# ----------------------------------

CXXLD=$(STATICTOOL) $(CXX)

bin_PROGRAMS=$(MEDITPROG)
EXTRA_DIST=picking.c 
EXTRA_PROGRAMS=ffmedit
ffmedit_SOURCES=animat.c cube.c image.c listnum.c mouse.c persp.c stream.c zaldy2.c bbfile.c dlists.c inmsh2.c \
	material.c normal.c tiles.c camera.c ellipse.c inout.c medit.c param.c psfile.c transform.c cenrad.c geometry.c \
	inout_morice.c menus.c parsar.c scene.c util.c clip.c gisfil.c inout_popenbinaire.c mesh.c parsop.c scissor.c vector.c \
	clipvol.c hash.c items.c mlists.c particle.c sftcpy.c view.c critip.c ilists.c keyboard.c morphing.c path.c status.c \
	zaldy1.c chrono.h extern.h image.h mesh.h eigenv.h grafic.h medit.h sproto.h


LDADD=picking.$(OBJEXT) ../libMesh/libMesh.a $(LIBSGLUT) -lm
AM_CPPFLAGS=-I$(srcdir)/../libMesh
BUILT_SOURCES=compil.date
BUILD_DATE="$(date --utc --date="@${SOURCE_DATE_EPOCH:-$(date +%s)}" +%Y-%m-%d)"
compil.date: $(ffmedit_SOURCES)
	echo "#define COMPIL " '"' $(BUILD_DATE) '(with ff++ $(VERSION))''"' > compil.date
#.PHONY: compil.date

# ALH - during a parallel make, we should make sure that picking.o is not used before being made by a different
# thread. If I understand things correctly, picking.c should not be compiled with the regular CFLAGS options.

ffmedit_DEPENDENCIES=picking.$(OBJEXT)
picking.$(OBJEXT):picking.c
	${CC} -c $< $(CNOFLAGS) -o $@ 

#pinking.no-optffmedit$(EXEEXT):picking.c
#	rm picking.$(OBJEXT)
#	$(MAKE) picking.o $(MEDITPROG) CFLAGS="$(CNOFLAGS)"
#	touch pinking.no-optffmedit$(EXEEXT)
#pinking.no-opt:picking.c
#	touch pinking.no-opt
#clean-local:
#	-rm pinking.no-opt pinking.no-optffmedit$(EXEEXT)

