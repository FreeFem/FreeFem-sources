BAMG=../../src/bamg/bamg
check:
	cd square;bamg=$(BAMG);export bamg;./dotest
	cd square;bamg=$(BAMG);export bamg;./dotest -q
	cd square;bamg=$(BAMG);export bamg;./doadapt
	cd test;bamg=$(BAMG);perl ./dotest.pl
