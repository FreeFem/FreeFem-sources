PARAM=-ne -nw -v 1000000

FF=../src/nw/FreeFem++
SRC=$(wildcard base/*.edp)

FFmpi=../src/mpi/ff-mpirun -np 4
SRCmpi=$(wildcard parallel/*.edp)

SRCplugins=$(wildcard plugins/*.edp)

.PHONY: all base parallel plugins

all: freefem++.pref base parallel plugins

freefem++.pref:
	echo loadpath += "../plugin/seq" > freefem++.pref

base:
	(for test in $(SRC); do \
		echo $$test; \
		$(FF) $$test $(PARAM) > $$test.log; \
	done)

parallel:
	(for test in $(SRCmpi); do \
		echo $$test; \
		$(FFmpi) $$test $(PARAM) > $$test.log; \
	done)

plugins:
	(for test in $(SRCplugins); do \
		echo $$test; \
		$(FF) $$test $(PARAM) > $$test.log; \
	done)

clean:
	rm base/*.log
	rm parallel/*.log
	rm plugins/*.log
