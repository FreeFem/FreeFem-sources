PARAM=-ne -nw -v 1000000
FF=../src/nw/FreeFem++
SRC=$(wildcard *.edp)

FFmpi=../src/mpi/ff-mpirun -np 4
SRCmpi=$(wildcard parallel/*.edp)

.PHONY: all base parallel


all: freefem++.pref base parallel

base:
	(for test in $(SRC); do \
		$(FF) $$test $(PARAM); \
	done)

freefem++.pref:
	echo loadpath += "../plugin/seq" > freefem++.pref

parallel:
	(for test in $(SRCmpi); do \
		$(FFmpi) $$test $(PARAM); \
	done)
