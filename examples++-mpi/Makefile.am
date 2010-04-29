# $Id$

TESTS=regtests.sh
TESTS_ENVIRONMENT=PROGLOC=$(top_srcdir)/src EXEEXT=$(EXEEXT) MPIPROG=$(MPIPROG)
EXTRA_DIST=*.edp *.idp  regtests.sh 

.cpp.$(DYLIB_SUFFIX): ff-c++
	 ../examples++-load//ff-c++ -auto  $< 
