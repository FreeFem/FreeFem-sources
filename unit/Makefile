FF=../src/nw/FreeFem++ -ne -nw
SRC=$(wildcard *.edp)

all:
	(for test in $(SRC); do \
		$(FF) $$test; \
	done)
