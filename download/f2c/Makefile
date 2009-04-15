
DIRPKG=../pkg
LISTINST=  ../bin/f2c ../bin/fort77 ../include/f2c.h ../lib/libf2c.a
DIRINSTALL=/usr/local
F2C_TAR=$(DIRPKG)/f2c.tar
all:compile

$(F2C_TAR):
	mkdir -p $(DIRPKG)
	$(WGET) # wget --passive-ftp ftp://netlib.bell-labs.com/netlib/f2c.tar
	curl http://netlib.sandia.gov/cgi-bin/netlib/netlibfiles.tar?filename=netlib/f2c -o "$(DIRPKG)/f2c.tar"
f2c:$(F2C_TAR)
	tar xvf $(F2C_TAR)
	gunzip -rf f2c/*
	gunzip -rf f2c/*/*
	cd f2c;unzip -d libf2c  libf2c.zip

f2c/libf2c/makefile: Makefile
	egrep -v 'ld -r|mv [$$]' <f2c/libf2c/makefile.u >f2c/libf2c/makefile

compile:f2c f2c/libf2c/makefile
	cd f2c/src; make -f makefile.u
	cd f2c/libf2c; make 

compile-10.4:f2c Makefile-MacOs
	cd f2c/src; make -f makefile.u
	cd f2c/libf2c; make -f ../../Makefile-MacOs OS=10.4
compile-10.5:f2c Makefile-MacOs
	cd f2c/src; make -f makefile.u
	cd f2c/libf2c; make -f ../../Makefile-MacOs OS=10.5


install:  ../lib ../include ../bin ../bin/fort77
	cp f2c/src/f2c ../bin
	cp f2c/src/f2c.h ../include
	cp f2c/libf2c/libf2c.a ../lib
	../bin/fort77 tt.f
	rm a.out

install-10.4: compile-10.4
	sudo make install-sudo
	make install
install-10.5: compile-10.5
	sudo make install-sudo
	make install
install-sudo: $(DIRINSTALL)/bin/fort77
	mkdir -p $(DIRINSTALL)/bin $(DIRINSTALL)/include $(DIRINSTALL)/lib
	cp f2c/src/f2c $(DIRINSTALL)/bin; 
	cp f2c/src/f2c.h $(DIRINSTALL)/include; 
	cp f2c/libf2c/libf2c.a $(DIRINSTALL)/lib 	 
	$(DIRINSTALL)/bin/fort77 tt.f	
	rm a.out



../bin/fort77: fort77.sed Makefile
	dd="`pwd`/.." ; \
	sed -e "s;@CC@;$(CC);"g \
	-e "s;@INC@;-I$$dd/include;g" \
         -e "s;@LLIBDIR@;-L$$dd/lib;g" \
         -e "s;@f2c@;$$dd/bin/f2c;g"  < fort77.sed >../bin/fort77
	chmod a+x ../bin/fort77

$(DIRINSTALL)/bin/fort77: fort77.sed Makefile
	dd="$(DIRINSTALL)/" ; \
	sed -e "s;@CC@;$(CC);"g \
	-e "s;@INC@;-I$$dd/include;g" \
         -e "s;@LLIBDIR@;-L$$dd/lib;g" \
         -e "s;@f2c@;$$dd/bin/f2c;g"  < fort77.sed >$(DIRINSTALL)/bin/fort77
	chmod a+x $(DIRINSTALL)/bin/fort77


../bin:	
	mkdir ../bin
../lib:	
	mkdir ../lib
../include:	
	mkdir ../include
clean:
	-rm -rf  f2c *~  a.out
veryclean: clean 
	-rm $(LISTINST) $(F2C_TAR)






