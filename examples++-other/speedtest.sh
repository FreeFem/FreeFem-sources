#!/bin/sh
# Speed benchmark for FreeFem++
# $Id$

# The FreeFem++ executable
ffexe=${PROGLOC}/nw/FreeFem++-nw${EXEEXT}
if test ! -x $ffexe
then
    echo $ffexe non existent
    exit 1
fi

# Write some build information into the trace file
echo ------------------------------------------- >> speedtest.out
date >> speedtest.out
echo $CXXFLAGS >> speedtest.out

# Run the actual test
$ffexe lap3-cpu.edp	\
		|grep -E --					\
		'-- lap (cholesky|CG|UMFPACK) +[0-9]+x[0-9]+'	\
		|tee -a speedtest.out
