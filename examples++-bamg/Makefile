BAMG=../../src/bamg/bamg
check:
	cd square;bamg=$(BAMG);export bamg;./dotest
	cd square;bamg=$(BAMG);export bamg;./dotest -q
	cd square;bamg=$(BAMG);export bamg;./doadapt
	cd test;bamg=$(BAMG);export bamg;perl ./dotest.pl
