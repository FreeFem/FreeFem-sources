include ../ff-flags

# working dirs
EXEDIR = #$()/executable

YAMS2DIR =$(abs_top_builddir)/download/yams/yams2$(yams_VERSION)
EXEDIR = $(YAMS2DIR)/executable
SRCDIR = $(YAMS2DIR)/sourceslib
OBJDIR = $(YAMS2DIR)/objects
ARCDIR = $(YAMS2DIR)/archives
DIRDIR = $(EXEDIR) $(OBJDIR) $(ARCDIR)
INCDIR = -I$(YAMS2DIR)/sourcesnew -I$(abs_top_builddir)/src/libMesh/
LDLDIR = -L$(abs_top_builddir)src/libMesh/ -lMesh 
VPATH  = $(SRCDIR)

# objects list
src    = $(wildcard $(SRCDIR)/*.c)
header = $(wildcard $(SRCDIR)/*.h)
objs   = $(patsubst $(SRCDIR)%,$(OBJDIR)%,$(src:.c=.o))
prog   = yams2
lib    = $(OBJDIR)/libyams2.a

#.SILENT:

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	$(CC) $(OPT64)  $(INCDIR) $(CFLAGS) -c $< -o $@

$(EXEDIR)/$(prog):$(DIRDIR) $(objs) 
	echo "#define COMPIL " '"' `date` '"' > $(SRCDIR)/compil.date
	#$(CC) -c $(CFLAGS) $(INCDIR) $(SRCDIR)/yams0.c -o $(OBJDIR)/yams0.o
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