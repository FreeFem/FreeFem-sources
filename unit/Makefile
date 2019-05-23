FF=../src/nw/FreeFem++ -ne -nw -v 1000000
SRC=$(wildcard *.edp)

all:
	(for test in $(SRC); do \
		$(FF) $$test; \
	done)
