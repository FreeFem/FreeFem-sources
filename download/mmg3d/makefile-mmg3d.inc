include ../ff-flags

# working dirs
MMG3DDIR = $(abs_top_builddir)/download/mmg3d/mmg3dlib
EXEDIR = $(MMG3DDIR)/bin/
SRCDIR = $(MMG3DDIR)/sources/
OBJDIR = $(MMG3DDIR)/objects/
ARCDIR = $(MMG3DDIR)/archives
DIRDIR = $(MMG3DDIR)/objects $(OBJDIR) $(ARCDIR)
INCDIR = -I$(MMG3DDIR)/sources/
LDLDIR = 
VPATH  = $(SRCDIR)

# objects list
src    = $(wildcard $(SRCDIR)/*.c)
header = $(wildcard $(SRCDIR)/*.h)
objs   = $(patsubst $(SRCDIR)%,$(OBJDIR)%,$(src:.c=.o))
lib    = $(OBJDIR)/libmmg3d.a


#.SILENT:

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(OPT64) $(INCDIR) $(CFLAGS) -c $< -o $@

$(lib): $(objs)
	@echo 'Fin Archivage' $@

$(objs):$(header)

$(DIRDIR):
	@[ -d $@ ] || mkdir $@

lib:	$(DIRDIR) $(objs)
	$(AR) $(ARFLAGS) $(lib) $(OBJDIR)/*.o  
	$(RANLIB) $(lib)
	@echo 'Fin Archivage' $@
				
clean:
	-rm $(objs) libmmg3d.a

tar: $(DIRDIR)
	tar czf $(ARCDIR)/$(prog).`date +"%Y.%m.%d"`.tgz sources makefile

target: $(EXEDIR)/$(prog)
