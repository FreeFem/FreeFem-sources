#!/bin/sh
# Testing procedure for MPI version of FreeFem++
# $Id$

NPROCS=2 REGTESTS=essai.edp ../regtests.sh
if test $? != 0
    then
    exit 1
fi
NPROCS=2 REGTESTS=schwarz.edp ../regtests.sh
if test $? != 0
    then
    exit 1
fi
NPROCS=2 REGTESTS=schwarz-b.edp ../regtests.sh
if test $? != 0
    then
    exit 1
fi
NPROCS=2 REGTESTS=schwarz-c.edp ../regtests.sh
if test $? != 0
    then
    exit 1
fi
