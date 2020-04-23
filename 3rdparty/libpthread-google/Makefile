include ff-flags
PACKAGE=../pkg/libpthread-google.tar.gz
SRCDIR= ./libpthread-google
libpthread-google: WHERE

ifeq ($(WIN32DLLTARGET),)
libpthread-google WHERE tag-install tag-tar $(PACKAGE):
	@echo "	Warning:not under windows do nothing, target: $@" 
else

libpthread-google: WHERE
	
WHERE:tag-install
		echo pthread-google  LD -L@DIR@/lib -lpthread-google  >../lib/WHERE.pthread-google
		echo pthread-google INCLUDE -I@DIR@/include/pthread-google >> ../lib/WHERE.pthread-google

tag-install:tag-tar
	cd $(SRCDIR)/src; \
	$(CC) -O -c *.c -I../include ;\
	$(AR) rcs ../../../lib/libpthread-google.a *.o 
	
	cp -r $(SRCDIR)/include/.  ../include/pthread-google
	touch $@

tag-tar: $(PACKAGE)
	tar zxvf  $(PACKAGE)
	touch $@
	
$(PACKAGE):
	../getall -o libpthread-google -a
	
endif
clean:
		-rm ff-flags tag-*
		-rm -rf libpthread-google*
		-rm FAIT $(FAIRE) 
		-rm ../lib/*pthread-google*
		-rm -rf ../include/pthread-google*


ff-flags: ../Makefile Makefile
	grep 'abs_top_builddir *=' ../Makefile > ff-flags
	grep 'CC *=' ../Makefile >> ff-flags
	grep 'CFLAGS *=' ../Makefile >> ff-flags
	grep 'LDFLAGS *=' ../Makefile >> ff-flags
	grep 'AR *=' ../Makefile >> ff-flags
	grep 'ARFLAGS *=' ../Makefile >> ff-flags
	grep 'RANLIB *=' ../Makefile >> ff-flags
	grep 'WIN32DLLTARGET *=' ../Makefile >> ff-flags
	grep 'WGET *=' ../Makefile >> ff-flags
	grep 'LIBS *=' ../Makefile >> ff-flags

