#!/bin/sh
# Testing procedure for MPI version of FreeFem++
# $Id$

# The MPI executable
ffexe=${PROGLOC}/mpi/FreeFem++-mpi${EXEEXT}

# Doing one parallel test
function dotest(){
    local cmd="mpirun -np 2 -machinefile machinefile $ffexe $1"
    echo $cmd
    eval $cmd
    if test $? != 0
	then
	exit 1
    fi
}

# Only if the executable was built
if test "$MPIPROG" != ""
then
    host=`hostname`
    echo $host>machinefile
    echo $host>>machinefile

    dotest essai.edp
    dotest schwarz.edp
    dotest schwarz-b.edp
    dotest schwarz-c.edp
fi
