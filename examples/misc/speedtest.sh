#!/bin/bash
# Speed benchmark for FreeFEM
# $Id$

# The FreeFEM executable
ffexe=${PROGLOC}/nw/FreeFem++-nw${EXEEXT}
if test ! -x $ffexe
then
    echo $ffexe non existent
    exit 1
fi

# Write some build information into the trace file
echo ------------------------------------------- >> speedtest.out
date >> speedtest.out
echo CXXFLAGS=$CXXFLAGS >> speedtest.out

# Run the actual test
$ffexe lap3-cpu.edp	\
		|grep -E --					\
		'-- lap (Cholesky|CG|UMFPACK|LU|Crout) +[0-9]+x[0-9]+'	\
		|tee -a speedtest.out
if test $PIPESTATUS != 0
    then
    exit 1
fi
