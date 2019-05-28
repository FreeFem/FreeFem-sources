FF=../src/nw/FreeFem++ -ne -nw -v 1000000
SRC=$(wildcard *.edp)

all: freefem++.pref
	(for test in $(SRC); do \
		$(FF) $$test; \
	done)

freefem++.pref:
	echo loadpath += "../plugin/seq" > freefem++.pref
