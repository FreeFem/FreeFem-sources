#!/bin/sh
# Copies Debian packages into a web server directory structure

# Any error stops the script
set -e

# $1 must point to the directory where the files will be copied
if test ! -d $1
    then
    echo "Usage: CopyToServer.sh <directory to copy files to>"
    exit 1
fi

# Directory structure we want to build
if test -d $1/dists
then
    rm -r $1/dists
fi
basedir='dists/packages/ff++/binary-i386'
mkdir -p $1/$basedir

# Server configuration file
cp apt-ftparchive.conf $1

# Copies all packages (potential improvement: we could spread
# architecture-independant packages into a separate directory).

cp ../../freefem++*.{dsc,deb,tar.gz} $1/$basedir

# Create package list
cd $1
apt-ftparchive packages ./$basedir > ./$basedir/Packages
gzip -c ./$basedir/Packages > ./$basedir/Packages.gz
