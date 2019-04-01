FF=FreeFem++ -ne
SRC=$(wildcard *.edp)

all:
	(for test in $(SRC); do \
		$(FF) $$test; \
	done)
