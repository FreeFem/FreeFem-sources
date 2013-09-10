include ../ff-flags

MSHMETDIR =  $(abs_top_builddir)/download/mshmet/mshmet$(mshmet_VERSION)
# working dirs
EXEDIR = #$(MSHMETDIR)/executable
SRCDIR = $(MSHMETDIR)/sourceslib
OBJDIR = $(MSHMETDIR)/objects/$(ARCHI)
ARCDIR = $(MSHMETDIR)/archives
DIRDIR = $(EXEDIR) $(OBJDIR) $(ARCDIR)
VPATH  = $(SRCDIR)
INCDIR = -I$(MSHMETDIR)/sourceslib/ -I$(abs_top_builddir)/src/libMesh/
LDLDIR = -L$(abs_top_builddir)/src/libMesh/ -lMesh

# objects list
src    = $(wildcard $(SRCDIR)/*.c)
header = $(wildcard $(SRCDIR)/*.h)
objs   = $(patsubst $(SRCDIR)%,$(OBJDIR)%,$(src:.c=.o))
prog   = mshmet
lib    = $(OBJDIR)/libmshmet.a

#.SILENT:

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(OPT64) $(INCDIR) $(CFLAGS) -c $< -o $@

$(EXEDIR)/$(prog):$(DIRDIR) $(objs)
	echo "#define COMPIL " '"' `date` '"' > $(SRCDIR)/compil.date
	#$(CC) -c $(CFLAGS) $(INCDIR) $(SRCDIR)/mshmet.c -o $(OBJDIR)/mshmet.o
	#$(CC) $(LDFLAGS) $(OPT64) $(LDLDIR) $(objs) -o $@ -lm
	$(AR) $(ARFLAGS) $(lib) $(objs)
	$(RANLIB) $(lib)

$(objs):$(header)

$(DIRDIR):
	@[ -d $@ ] || mkdir $@

clean:
	-rm $(objs) $(EXEDIR)/$(prog)

tar:$(DIRDIR)
	tar czf $(ARCDIR)/$(prog).`date +"%Y.%m.%d"`.tgz sources makefile

target: $(EXEDIR)/$(prog)
