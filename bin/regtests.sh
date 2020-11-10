#!/usr/bin/env bash
# Runs all regression tests on all compiled FreeFem++ versions
# ------------------------------------------------------------
MPIRUN=`awk  '$1 =="MPIRUN" {print $3}' Makefile`
# Antoine Le Hyaric - LJLL Paris 6 - lehyaric@ann.jussieu.fr - 19/5/04
# $Id$

# To run one set of tests on one executable
# $1=program to run
# $2=tag for trace file
# $3=EDP script to run
function dotest(){

    # Running FreeFem++ on regtests.edp (specific to regression
    # tests), otherwise on all.edp.
    echo regtests.sh: running $1 $3, result in regtests-$2.log
    $1 $3|tee regtests-$2.log
    if test $PIPESTATUS != 0
	then
	exit 1
    fi
}

# For the example++-load tests
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:."

# In visual checks, we can run even the most invasive programs
script=$REGEDP
if test "$VISUALCHECK" = yes -a "$REGEDP" = regtests.edp
then
    script=all.edp
fi

# Number of processors in parallel mode
if test "$NPROCS" != ""
    then
    nprocs=$NPROCS
else
    nprocs=1
fi

# Do not test windowed programs by default, because their windows are
# too invasive.

if test "$VISUALCHECK" = yes
then
    export PATH="${PROGLOC}/nw/:$PATH";  dotest FreeFem++${EXEEXT} std $script
fi

if test $nprocs = 1
    then
    dotest ${PROGLOC}/nw/FreeFem++-nw${EXEEXT} nw $script
fi

if test "${X11PROG}" != "" -a "${VISUALCHECK}" = "yes"
then
    dotest ${PROGLOC}/x11/FreeFem++-x11${EXEEXT} x11 $script
fi

if test "${GLXPROG}" != "" -a "${VISUALCHECK}" = "yes"
then
    dotest ${PROGLOC}/glx/FreeFem++-glx${EXEEXT} glx $script
fi

if test "${AGLPROG}" != "" -a "${VISUALCHECK}" = "yes"
then
    dotest ${PROGLOC}/agl/FreeFem++-agl${EXEEXT} agl $script
fi

if test "${MPIPROG}" != ""
    then
    mpich=`${MPIRUN=mpirun} -h 2>&1 |grep mpich  |wc  -l`
    host=`hostname`
    echo $host>machinefile
    echo $host>>machinefile
    if [ $mpich -ne 0 ] ; then
    dotest "${MPIRUN} -np $nprocs -machinefile machinefile ${PROGLOC}/mpi/FreeFem++-mpi${EXEEXT}" mpi $script
    else
	[[ -f "$(which lamboot 2>/dev/null)" ]] && lamboot	
     dotest "${MPIRUN} -np $nprocs ${PROGLOC}/mpi/FreeFem++-mpi${EXEEXT}" mpi $script
    fi
fi

if test "${IDEPROG}" != "" -a "${VISUALCHECK}" = "yes"
then
    dotest ${PROGLOC}/ide/FreeFem++-cs${EXEEXT} ide $script
fi
