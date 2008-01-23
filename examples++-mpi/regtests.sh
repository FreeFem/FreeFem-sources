#!/bin/sh
# Testing procedure for MPI version of FreeFem++
# $Id$

NPROCS=2 REGEDP=essai.edp ../regtests.sh
if test $? != 0
    then
    exit 1
fi
NPROCS=2 REGEDP=schwarz.edp ../regtests.sh
if test $? != 0
    then
    exit 1
fi
NPROCS=2 REGEDP=schwarz-b.edp ../regtests.sh
if test $? != 0
    then
    exit 1
fi
NPROCS=2 REGEDP=schwarz-c.edp ../regtests.sh
if test $? != 0
    then
    exit 1
fi
NPROCS=4 REGEDP=mortar-DN-4-mpi.edp ../regtests.sh
if test $? != 0
    then
    exit 1
fi

