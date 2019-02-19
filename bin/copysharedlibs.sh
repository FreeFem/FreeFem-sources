#!/bin/sh
# Find out which shared libs an executable needs and copy them
# Antoine Le Hyaric - LJLL Paris 6 - lehyaric@ann.jussieu.fr - 22/11/04
# $Id$

# $1=executable to analyze
if test ! -x $1
then
    echo $1 is not an executable
    exit 1
fi

# $2= where to copy shared libs
if test ! -d $2
then
    echo $2 is not a directory
fi

# List all shared libs
libs=`ldd $1|awk '{print $3}'`
if test "$libs" != "dynamic" -a "$libs" != ""
then
    cp $libs $2
fi
